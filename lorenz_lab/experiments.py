"""Воспроизводимые эксперименты. Числа отчёта берутся только отсюда."""
import csv
import json
import platform
from pathlib import Path
from time import perf_counter
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from lorenz import Parameters, rhs, jacobian, equilibria
from numerical_methods import integrate, time_grid

DEFAULT_CONFIG = Path(__file__).with_name('config.json')

def load_config(path=DEFAULT_CONFIG):
    c = json.loads(Path(path).read_text(encoding='utf-8'))
    Parameters(c['sigma'], c['rho'], c['beta'])
    for key in ['initial_state', 'perturbation']:
        a = np.asarray(c[key], dtype=float)
        if a.shape != (3,) or not np.all(np.isfinite(a)):
            raise ValueError(key+' должен содержать три конечных числа')
    if not c['t_start'] < c['short_end'] <= c['t_end']:
        raise ValueError('Короткий интервал должен лежать внутри основного')
    if not c['t_start'] <= c['transient_end'] < c['t_end']:
        raise ValueError('Граница переходного процесса должна лежать внутри интервала')
    if np.linalg.norm(np.array(c['perturbation'])) == 0:
        raise ValueError('Для чувствительности требуется ненулевое возмущение')
    time_grid(c['t_start'], c['t_end'], c['step'])
    for h in c['steps']+c['euler_refinement_steps']: time_grid(c['t_start'], c['short_end'], h)
    if not c['t_start'] <= c['slope_window'][0] < c['slope_window'][1] <= c['t_end']:
        raise ValueError('Окно оценки наклона должно лежать в интервале расчёта')
    if c['timing_repeats'] < 1 or c['animation_frames'] < 2 or c['animation_fps'] < 1 or c['dpi'] < 150:
        raise ValueError('Некорректные настройки повторений, GIF или dpi')
    for k in ['rk45_rtol','rk45_atol','reference_rtol','reference_atol']:
        if not np.isfinite(c[k]) or c[k] <= 0: raise ValueError('Допуски должны быть положительными')
    return c

def adaptive(fun, initial, start, end, rtol, atol, method='RK45'):
    # t_eval не передаём: sol.t — внутренняя сетка принятых шагов.
    sol = solve_ivp(fun, (start, end), initial, method=method,
                    rtol=rtol, atol=atol, dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

def timed(call, repeats):
    durations=[]; result=None
    for _ in range(repeats):
        t0=perf_counter(); result=call(); durations.append(perf_counter()-t0)
    return result, float(np.median(durations))

def write_csv(path, rows):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def run_experiments(c, output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    p=Parameters(c['sigma'],c['rho'],c['beta']); fun=lambda t,u:rhs(t,u,p)
    t0,t1,h=c['t_start'],c['t_end'],c['step'];u0=np.array(c['initial_state'])
    grid=time_grid(t0,t1,h); reference=adaptive(fun,u0,t0,t1,c['reference_rtol'],c['reference_atol'])
    ref=reference.sol(grid).T
    solutions={};metrics=[]
    for method in ['euler','rk4','RK45']:
        if method=='RK45':
            def call():
                sol=adaptive(fun,u0,t0,t1,c['rk45_rtol'],c['rk45_atol'])
                return sol,sol.sol(grid).T
            (sol,values),seconds=timed(call,c['timing_repeats']);steps=len(sol.t)-1;nfev=sol.nfev
        else:
            sol,seconds=timed(lambda:integrate(fun,u0,t0,t1,h,method),c['timing_repeats'])
            values=sol.y;steps=sol.n_steps;nfev=sol.nfev
        error=np.linalg.norm(values-ref,axis=1);solutions[method]=values
        short=grid<=c['short_end']
        metrics.append(dict(method=method,max_error=float(max(error)),final_error=float(error[-1]),
                            short_max_error=float(max(error[short])),seconds_median=seconds,
                            n_steps=steps,nfev=nfev,n_output=len(grid)))
    # Независимая проверка эталона только на коротком интервале.
    short_ref=adaptive(fun,u0,t0,c['short_end'],c['reference_rtol'],c['reference_atol'])
    # Контрольный rtol выше нижнего предела SciPy 100*eps(float64).
    cross=adaptive(fun,u0,t0,c['short_end'],3e-14,1e-15,'DOP853')
    tight=adaptive(fun,u0,t0,c['short_end'],3e-14,1e-15)
    short_grid=time_grid(t0,c['short_end'],min(c['steps']))
    ref_short=short_ref.sol(short_grid).T
    reference_cross_error=float(np.max(np.linalg.norm(ref_short-cross.sol(short_grid).T,axis=1)))
    reference_tight_error=float(np.max(np.linalg.norm(ref_short-tight.sol(short_grid).T,axis=1)))
    convergence=[]
    for method in ['euler','rk4']:
        previous=None
        for step in sorted(c['steps'],reverse=True):
            try:
                sol=integrate(fun,u0,t0,c['short_end'],step,method)
                err=np.linalg.norm(sol.y-short_ref.sol(sol.t).T,axis=1)
                row=dict(method=method,h=step,max_error=float(max(err)),final_error=float(err[-1]),
                         observed_order='',status='ok')
                if previous:
                    row['observed_order']=float(np.log(previous['max_error']/row['max_error'])/np.log(previous['h']/step))
                previous=row
            except FloatingPointError:
                row=dict(method=method,h=step,max_error='',final_error='',observed_order='',status='nonfinite')
                previous=None
            convergence.append(row)
    refinement=[]
    for step in c['euler_refinement_steps']:
        sol=integrate(fun,u0,t0,c['short_end'],step,'euler')
        err=np.linalg.norm(sol.y-short_ref.sol(sol.t).T,axis=1)
        refinement.append(dict(h=step,max_error=float(max(err))))
    refinement_order=float(np.polyfit(np.log([r['h'] for r in refinement[-3:]]),np.log([r['max_error'] for r in refinement[-3:]]),1)[0])
    write_csv(out/'euler_refinement.csv',refinement)
    perturbed=integrate(fun,u0+np.array(c['perturbation']),t0,t1,h,'rk4').y
    distance=np.linalg.norm(solutions['rk4']-perturbed,axis=1)
    mask=(grid>=c['slope_window'][0])&(grid<=c['slope_window'][1])&(distance>0)
    slope,intercept=np.polyfit(grid[mask],np.log(distance[mask]),1)
    log_pred=slope*grid[mask]+intercept
    denom=np.sum((np.log(distance[mask])-np.mean(np.log(distance[mask])))**2)
    r2=float(1-np.sum((np.log(distance[mask])-log_pred)**2)/denom) if denom>0 else None
    # Повтор чувствительности с половинным шагом, сравниваем только оценку и расстояния.
    fine=integrate(fun,u0,t0,t1,h/2,'rk4')
    fine_p=integrate(fun,u0+np.array(c['perturbation']),t0,t1,h/2,'rk4')
    fd=np.linalg.norm(fine.y-fine_p.y,axis=1)
    fm=(fine.t>=c['slope_window'][0])&(fine.t<=c['slope_window'][1])&(fd>0)
    fine_slope=float(np.polyfit(fine.t[fm],np.log(fd[fm]),1)[0])
    eqrows=[]
    for name,state in equilibria(p).items():
        eigen=np.linalg.eigvals(jacobian(state,p));reals=np.real(eigen)
        # 1e-10 — фильтр почти нулевых Re(lambda), не физический порог.
        stability='неустойчиво' if np.any(reals>1e-10) else 'асимптотически устойчиво' if np.all(reals< -1e-10) else 'линеаризация не решает вопрос'
        eqrows.append(dict(name=name,state=state.tolist(),eigenvalues=[{'real':float(z.real),'imag':float(z.imag)} for z in eigen],stability=stability))
    # Результаты и промежуточные массивы сохраняются для независимого анализа.
    write_csv(out/'metrics.csv',metrics);write_csv(out/'convergence.csv',convergence)
    errors={k:np.linalg.norm(v-ref,axis=1) for k,v in solutions.items()}
    np.savez_compressed(out/'trajectories.npz',t=grid,reference=ref,euler=solutions['euler'],rk4=solutions['rk4'],rk45=solutions['RK45'],perturbed=perturbed,distance=distance,fine_t=fine.t,fine_distance=fd)
    orders={}
    for method in ['euler','rk4']:
        valid=[r for r in convergence if r['method']==method and r['status']=='ok']
        chosen=valid[-3:]
        orders[method]=float(np.polyfit(np.log([r['h'] for r in chosen]),np.log([r['max_error'] for r in chosen]),1)[0]) if len(chosen)>=2 else None
    summary=dict(config=c,versions={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
                 metrics=metrics,convergence=convergence,euler_refinement=refinement,euler_refinement_order=refinement_order,orders_last_three=orders,equilibria=eqrows,
                 reference=dict(rtol=c['reference_rtol'],atol=c['reference_atol'],n_steps=len(reference.t)-1,nfev=reference.nfev,
                                short_dop853_disagreement=reference_cross_error,short_tighter_rk45_disagreement=reference_tight_error),
                 sensitivity=dict(slope=float(slope),intercept=float(intercept),r_squared=r2,window=c['slope_window'],
                                  fine_slope=fine_slope,initial_distance=float(distance[0]),final_distance=float(distance[-1]),
                                  first_distance_above_one=float(grid[np.flatnonzero(distance>=1)[0]]) if np.any(distance>=1) else None),
                 all_standard_trajectories_finite=bool(all(np.all(np.isfinite(x)) for x in [*solutions.values(),ref,perturbed])))
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    return dict(summary=summary,t=grid,solutions=solutions,reference=ref,errors=errors,
                convergence=convergence,perturbed=perturbed,distance=distance,equilibria=eqrows)
