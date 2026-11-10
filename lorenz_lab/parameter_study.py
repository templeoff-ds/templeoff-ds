"""Дополнительные эксперименты для защиты выбора численных констант.

Запуск: python parameter_study.py. Исходная модель и её результаты не меняются.
Диапазоны альтернатив явно заданы в constants-study.json, а не спрятаны в коде.
"""
import ast
import json
import platform
from pathlib import Path
from time import perf_counter
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lorenz import Parameters, rhs, jacobian, equilibria
from numerical_methods import integrate, time_grid
from experiments import load_config, adaptive, write_csv

ROOT=Path(__file__).resolve().parent

def slope_result(t,d,window):
    mask=(t>=window[0])&(t<=window[1])&(d>0)
    if np.count_nonzero(mask)<2:
        return dict(slope=None,r_squared=None,d_min=float(d.min()),d_max=float(d.max()))
    a,b=np.polyfit(t[mask],np.log(d[mask]),1)
    residual=np.sum((np.log(d[mask])-(a*t[mask]+b))**2)
    total=np.sum((np.log(d[mask])-np.log(d[mask]).mean())**2)
    return dict(slope=float(a),r_squared=float(1-residual/total) if total>0 else None,
                d_min=float(d[mask].min()),d_max=float(d[mask].max()))

def study():
    c=load_config();q=json.loads((ROOT/'constants-study.json').read_text())
    out=ROOT/'results'/'constants_study';out.mkdir(parents=True,exist_ok=True)
    p=Parameters(c['sigma'],c['rho'],c['beta']);u=np.array(c['initial_state']);fun=lambda t,x:rhs(t,x,p)
    t=time_grid(c['t_start'],c['t_end'],c['step'])
    base=integrate(fun,u,c['t_start'],c['t_end'],c['step'])
    pert=integrate(fun,u+np.array(c['perturbation']),c['t_start'],c['t_end'],c['step'])
    d=np.linalg.norm(pert.y-base.y,axis=1)
    ref=adaptive(fun,u,c['t_start'],c['t_end'],c['reference_rtol'],c['reference_atol'])
    result={'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'baseline':c,'study_config':q}
    def save(name,rows):
        result[name]=rows;write_csv(out/(name+'.csv'),rows)
    rows=[]
    for key,values in q['model_variants'].items():
        for value in values:
            pp=Parameters(**{**dict(sigma=p.sigma,rho=p.rho,beta=p.beta),key:value})
            ss=adaptive(lambda t,x:rhs(t,x,pp),u,c['t_start'],c['t_end'],c['reference_rtol'],c['reference_atol'])
            ys=ss.sol(t).T;tail=ys[t>=c['transient_end']]
            rows.append(dict(parameter=key,value=value,max_state_norm=float(np.max(np.linalg.norm(ys,axis=1))),
                             tail_std_x=float(np.std(tail[:,0])),
                             equilibria_count=len(equilibria(pp)),max_real_eigen_Cplus=float(np.max(np.real(np.linalg.eigvals(jacobian(equilibria(pp)['C+'],pp))))) if pp.rho>1 else None))
    save('model_parameters',rows)
    rows=[]
    for initial in q['initial_variants']:
        sol=adaptive(fun,initial,c['t_start'],c['t_end'],c['reference_rtol'],c['reference_atol'])
        ys=sol.sol(t).T
        rows.append(dict(initial=str(initial),initial_speed=float(np.linalg.norm(fun(t[0],initial))),max_norm=float(np.max(np.linalg.norm(ys,axis=1))),final_norm=float(np.linalg.norm(ys[-1]))))
    save('initial_states',rows)
    # Контрольная сетка одинакова для всех шагов: мелкая сетка исходного исследования.
    rows=[]
    for end in q['short_horizons']:
        errors=[]
        for h in sorted(c['steps'],reverse=True):
            sol=integrate(fun,u,c['t_start'],end,h,'rk4') if np.isclose(round((end-c['t_start'])/h)*h,end-c['t_start']) else None
            if sol is None:continue
            error=float(np.max(np.linalg.norm(sol.y-ref.sol(sol.t).T,axis=1)))
            errors.append((h,error))
        chosen=errors[-3:]
        order=float(np.polyfit(np.log([r[0] for r in chosen]),np.log([r[1] for r in chosen]),1)[0])
        rows.append(dict(end=end,order_last_three=order,error_h001=next(e for h,e in errors if h==c['step'])))
    save('short_intervals',rows)
    horizons=[]
    for end in q['horizons']:
        if end <= c['t_end']:
            error=float(np.max(np.linalg.norm(base.y[t<=end]-ref.sol(t[t<=end]).T,axis=1)))
            distance=float(d[np.searchsorted(t,end)])
        else:
            longer=integrate(fun,u,c['t_start'],end,c['step'])
            longer_p=integrate(fun,u+np.array(c['perturbation']),c['t_start'],end,c['step'])
            longer_ref=adaptive(fun,u,c['t_start'],end,c['reference_rtol'],c['reference_atol'])
            error=float(np.max(np.linalg.norm(longer.y-longer_ref.sol(longer.t).T,axis=1)))
            distance=float(np.linalg.norm(longer.y[-1]-longer_p.y[-1]))
        n=int(round((end-c['t_start'])/c['step']))
        horizons.append(dict(end=end,max_rk4_error=error,distance_at_end=distance,n_steps=n,nfev=4*n))
    save('horizons',horizons)
    save('distance_thresholds',[dict(threshold=threshold,first_crossing=float(t[np.flatnonzero(d>=threshold)[0]]) if np.any(d>=threshold) else None) for threshold in q['distance_thresholds']])
    short_grid=time_grid(c['t_start'],c['short_end'],min(c['steps']))
    control=adaptive(fun,u,c['t_start'],c['short_end'],q['control_rtol'],q['control_atol'],'DOP853').sol(short_grid).T
    rows=[]
    for rtol,atol in q['tolerance_pairs']:
        seconds=[]
        for _ in range(c['timing_repeats']):
            begin=perf_counter();sol=adaptive(fun,u,c['t_start'],c['short_end'],rtol,atol)
            ys=sol.sol(short_grid).T;seconds.append(perf_counter()-begin)
        rows.append(dict(rtol=rtol,atol=atol,max_error=float(np.max(np.linalg.norm(ys-control,axis=1))),nfev=sol.nfev,n_steps=len(sol.t)-1,seconds_median=float(np.median(seconds))))
    save('tolerances',rows)
    rows=[]
    for delta in q['perturbations']:
        perturbed_u=u.copy();perturbed_u[0]+=delta
        ys=integrate(fun,perturbed_u,c['t_start'],c['t_end'],c['step']).y
        dist=np.linalg.norm(ys-base.y,axis=1);reached=np.flatnonzero(dist>=1)
        rows.append(dict(requested_delta=delta,actual_delta=float(perturbed_u[0]-u[0]),first_above_one=float(t[reached[0]]) if reached.size else None,**slope_result(t,dist,c['slope_window'])))
    save('perturbations',rows)
    save('windows',[dict(start=a,end=b,**slope_result(t,d,[a,b])) for a,b in q['slope_windows']])
    save('discard_times',[dict(discard=begin,retained_fraction=float(np.mean(t>=begin)),retained_points=int(np.count_nonzero(t>=begin)),std_x=float(np.std(base.y[t>=begin,0]))) for begin in q['discard_times']])
    # Эффект числа последних точек, а не подбор лучшего наклона.
    existing=json.loads((ROOT/'results'/'summary.json').read_text())
    if existing['config'] != c:
        raise ValueError('Конфигурация изменилась: сначала повторите python main.py')
    rows=[]
    for method in ['euler','rk4']:
        good=[r for r in existing['convergence'] if r['method']==method and r['status']=='ok']
        for count in q['fit_points']:
            if count>len(good):continue
            chosen=good[-count:]
            rows.append(dict(method=method,point_count=count,h_max=max(r['h'] for r in chosen),h_min=min(r['h'] for r in chosen),order=float(np.polyfit(np.log([r['h'] for r in chosen]),np.log([r['max_error'] for r in chosen]),1)[0])))
    save('fit_points',rows)
    rows=[]
    for repeats in q['timing_repeat_counts']:
        durations=[]
        for _ in range(repeats):
            begin=perf_counter();integrate(fun,u,c['t_start'],c['t_end'],c['step']);durations.append(perf_counter()-begin)
        rows.append(dict(repeats=repeats,median_seconds=float(np.median(durations)),min_seconds=float(min(durations)),max_seconds=float(max(durations))))
    save('timings',rows)
    state=np.array([2.3,-1.4,8.1]);rows=[]
    for eps in q['jacobian_eps']:
        fd=np.column_stack([(fun(0,state+eps*e)-fun(0,state-eps*e))/(2*eps) for e in np.eye(3)])
        rows.append(dict(eps=eps,max_entry_error=float(np.max(np.abs(fd-jacobian(state,p))))))
    save('jacobian_difference',rows)
    rows=[]
    for threshold in q['eigen_thresholds']:
        for name,state in equilibria(p).items():
            es=np.real(np.linalg.eigvals(jacobian(state,p)))
            label='unstable' if np.any(es>threshold) else 'stable' if np.all(es < -threshold) else 'inconclusive'
            rows.append(dict(threshold=threshold,equilibrium=name,classification=label))
    save('eigen_thresholds',rows)
    rows=[]
    for rtol,atol in q['grid_tolerances']:
        # Первый случай имеет только двоичное округление, второй — намеренное несовпадение.
        for label,start,end,h in [('binary_roundoff',0,.3,.1),('near_noninteger',0,1+1e-12,.01),('noninteger',0,1+1e-7,.01)]:
            n=round((end-start)/h)
            rows.append(dict(rtol=rtol,atol=atol,case=label,length_disagreement=float(abs(n*h-(end-start))),accepted=bool(np.isclose(n*h,end-start,rtol=rtol,atol=atol))))
    save('grid_tolerances',rows)
    # Ошибки, с которыми сравниваются численные пороги автоматических тестов.
    rows=[]
    exp_error=abs(integrate(lambda t,x:x,[1],0,1,.01).y[-1,0]-np.e)
    nonauto=abs(integrate(lambda t,x:t*x,[1],0,1,.01).y[-1,0]-np.exp(.5))
    check=integrate(fun,u,0,2,.001);control_dense=adaptive(fun,u,0,2,q['control_rtol'],q['control_atol'],'DOP853')
    short_error=float(np.max(np.linalg.norm(check.y-control_dense.sol(check.t).T,axis=1)))
    for name,error,threshold in [('rk4_exp',exp_error,3e-10),('rk4_nonautonomous',nonauto,1e-8),('rk4_lorenz',short_error,2e-6)]:
        rows.append(dict(test=name,measured_error=float(error),threshold=threshold,threshold_over_error=float(threshold/error)))
    save('test_error_margins',rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    rows=[]
    for dpi in q['dpi_variants']:
        fig,ax=plt.subplots(figsize=(6.4,4.8));ax.plot(t,base.y[:,0]);ax.set(xlabel='t',ylabel='x',title=f'Контроль вывода: dpi={dpi}')
        begin=perf_counter();filename=out/f'dpi_{dpi}.png';fig.savefig(filename,dpi=dpi);elapsed=perf_counter()-begin;plt.close(fig)
        from PIL import Image
        with Image.open(filename) as im:width,height=im.size
        rows.append(dict(dpi=dpi,width_px=width,height_px=height,bytes=filename.stat().st_size,save_seconds=elapsed))
    save('image_resolution',rows)
    rows=[]
    for nframes in q['frame_variants']:
        idx=np.unique(np.linspace(1,len(t),nframes).astype(int))-1
        for fps in q['fps_variants']:
            rows.append(dict(frames=len(idx),fps=fps,duration_seconds=float(len(idx)/fps),max_simulation_time_gap=float(np.max(np.diff(t[idx]))),all_integrator_points_still_in_line=True))
    save('animation_sampling',rows)
    # Сводные рисунки не затрагивают исходные графики.
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    r=result['tolerances'];axes[0,0].loglog([x['nfev'] for x in r],[x['max_error'] for x in r],'o-');axes[0,0].set(xlabel='Вызовы правой части',ylabel='Ошибка на [0,2]',title='Допуски: точность и стоимость')
    r=result['perturbations'];axes[0,1].semilogx([x['actual_delta'] for x in r if x['actual_delta']>0],[x['first_above_one'] if x['first_above_one'] is not None else np.nan for x in r if x['actual_delta']>0],'o-');axes[0,1].set(xlabel='Фактическое возмущение',ylabel='Первый t с d≥1',title='Меньшая δ откладывает расхождение')
    r=result['windows'];axes[1,0].bar(range(len(r)),[x['slope'] for x in r]);axes[1,0].set_xticks(range(len(r)),[f"{x['start']}–{x['end']}" for x in r]);axes[1,0].set(xlabel='Окно регрессии',ylabel='Наклон',title='Зависимость оценки от окна')
    r=result['jacobian_difference'];axes[1,1].loglog([x['eps'] for x in r],[max(x['max_entry_error'],np.finfo(float).tiny) for x in r],'o-');axes[1,1].set(xlabel='Приращение ε',ylabel='max ошибки J',title='Разностный Якобиан: слишком малое ε вредно')
    for ax in axes.flat:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'constant_tradeoffs.png',dpi=c['dpi']);plt.close(fig)
    (out/'study.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Результаты проверки констант:',out)
    return result

if __name__=='__main__':study()
