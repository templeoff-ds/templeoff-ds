"""Графики и GIF для запуска на сервере без графического интерфейса."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from pathlib import Path

COLORS={'euler':'#dc6848','rk4':'#246cab','RK45':'#3e9a68'}

def setup_style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.grid':True,'grid.alpha':.23,'lines.linewidth':1.1,
                         'axes.spines.top':False,'axes.spines.right':False})

def render_figures(data,c,output):
    setup_style();out=Path(output);out.mkdir(parents=True,exist_ok=True)
    t=data['t'];y=data['solutions']['rk4'];yp=data['perturbed'];keep=t>=c['transient_end']
    tag=f"σ={c['sigma']:g}; ρ={c['rho']:g}; β={c['beta']:.6g}; h={c['step']:g}"
    def save(fig,name):
        fig.savefig(out/name,dpi=c['dpi'],bbox_inches='tight');plt.close(fig)
    fig=plt.figure(figsize=(9,7));ax=fig.add_subplot(projection='3d')
    ax.plot(*y[keep].T,color=COLORS['rk4'],lw=.55,label='RK4 после переходного процесса')
    ax.set(xlabel='x',ylabel='y',zlabel='z',title='Численная траектория на аттракторе Лоренца\n'+tag);ax.legend();save(fig,'lorenz_attractor.png')
    fig,axes=plt.subplots(3,1,figsize=(10,7),sharex=True)
    for i,ax in enumerate(axes):ax.plot(t,y[:,i],color=['#246cab','#3e9a68','#dc6848'][i]);ax.set_ylabel('xyz'[i])
    axes[0].set_title('Временные зависимости, RK4\n'+tag);axes[-1].set_xlabel('t');fig.tight_layout();save(fig,'time_series.png')
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,(i,j) in zip(axes,[(0,1),(0,2),(1,2)]):
        ax.plot(y[keep,i],y[keep,j],lw=.5);ax.set(xlabel='xyz'[i],ylabel='xyz'[j])
    fig.suptitle('Фазовые портреты после переходного процесса\n'+tag);fig.tight_layout();save(fig,'phase_portraits.png')
    fig,axes=plt.subplots(3,1,figsize=(11,10))
    for method,ys in data['solutions'].items():
        axes[0].plot(t,ys[:,0],label=method,color=COLORS[method]);axes[1].plot(t,ys[:,0],label=method,color=COLORS[method])
        axes[2].semilogy(t,np.maximum(data['errors'][method],1e-16),label=method,color=COLORS[method])
    for ax in axes:ax.set_xlabel('t');ax.legend()
    axes[0].set(xlim=(c['t_start'],c['short_end']),ylabel='x',title='Сравнение на коротком интервале')
    axes[1].set(ylabel='x',title='Сравнение на полном интервале: поточечное совпадение теряется')
    axes[2].set(ylabel='‖u − u_ref‖₂',title='Ошибка относительно численного эталона RK45 (не точное решение)')
    fig.suptitle(tag);fig.tight_layout();save(fig,'methods_comparison.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for ax,method,order in zip(axes,['euler','rk4'],[1,4]):
        rows=[r for r in data['convergence'] if r['method']==method and r['status']=='ok']
        h=np.array([r['h'] for r in rows]);err=np.array([r['max_error'] for r in rows])
        ax.loglog(h,err,'o-',label='Измеренный максимум ошибки',color=COLORS[method])
        if method=='euler':
            extra=data['summary']['euler_refinement']
            ax.loglog([r['h'] for r in extra],[r['max_error'] for r in extra],'s:',label='Дополнительное измельчение')
        ax.loglog(h,err[-1]*(h/h[-1])**order,'--',label=f'Ориентир h^{order}')
        ax.set(xlabel='Шаг h',ylabel='max ‖u − u_ref‖₂',title=method+f", p≈{data['summary']['orders_last_three'][method]:.3f} (три малых шага)");ax.legend(fontsize=8)
    fig.suptitle(f"Сходимость на [{c['t_start']:g}, {c['short_end']:g}]; "+tag.rsplit('; h=',1)[0]);fig.tight_layout();save(fig,'step_dependence.png')
    fig=plt.figure(figsize=(12,10));ax=fig.add_subplot(221,projection='3d')
    ax.plot(*y.T,lw=.5,label='u₀');ax.plot(*yp.T,lw=.5,label='u₀+δ');ax.set(xlabel='x',ylabel='y',zlabel='z',title='Две близкие траектории');ax.legend()
    ax=fig.add_subplot(222);ax.plot(t,data['distance']);ax.set(xlabel='t',ylabel='d(t)',title='Евклидово расстояние')
    ax=fig.add_subplot(223);ax.semilogy(t,np.maximum(data['distance'],1e-17),label='h');
    with np.load(out/'trajectories.npz') as stored:ax.semilogy(stored['fine_t'],np.maximum(stored['fine_distance'],1e-17),label='h/2',alpha=.7)
    ax.set(xlabel='t',ylabel='d(t), log-шкала',title='Повтор с половинным шагом');ax.legend()
    ax=fig.add_subplot(224);ss=data['summary']['sensitivity'];mask=(t>=ss['window'][0])&(t<=ss['window'][1])
    ax.plot(t[mask],np.log(data['distance'][mask]),label='ln d(t)');ax.plot(t[mask],ss['intercept']+ss['slope']*t[mask],'--',label=f"Регрессия, a={ss['slope']:.3f}")
    ax.set(xlabel='t',ylabel='ln d(t)',title='Наклон на выбранном окне, не точный показатель Ляпунова');ax.legend(fontsize=8)
    fig.suptitle('Чувствительность к начальному состоянию\n'+tag);fig.tight_layout();save(fig,'sensitivity.png')
    fig=plt.figure(figsize=(11,5));ax=fig.add_subplot(121,projection='3d');ax.plot(*y[keep].T,alpha=.45,lw=.4)
    ax2=fig.add_subplot(122)
    for row in data['equilibria']:
        ax.scatter(*row['state'],s=50,label=row['name']);es=[complex(v['real'],v['imag']) for v in row['eigenvalues']]
        ax2.scatter(np.real(es),np.imag(es),s=45,label=row['name'],marker='x' if row['name']=='C-' else 'o')
    ax.set(xlabel='x',ylabel='y',zlabel='z',title='Равновесия и траектория');ax.legend();ax2.axvline(0,c='k',lw=.8)
    ax2.set(xlabel='Re λ',ylabel='Im λ',title='Спектры Якобиана');ax2.legend();fig.suptitle(tag);fig.tight_layout();save(fig,'equilibria.png')

def render_animation(data,c,output):
    setup_style();out=Path(output)
    t=data['t'];y=data['solutions']['rk4'];yp=data['perturbed']
    fig=plt.figure(figsize=(6.4,5.2));ax=fig.add_subplot(projection='3d')
    lo=np.minimum(y.min(axis=0),yp.min(axis=0))-2;hi=np.maximum(y.max(axis=0),yp.max(axis=0))+2
    ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]),xlabel='x',ylabel='y',zlabel='z')
    ax.set_title('Лоренц: две близкие траектории, RK4\n'+f"σ={c['sigma']:g}, ρ={c['rho']:g}, β={c['beta']:.4g}, h={c['step']:g}")
    line1,=ax.plot([],[],[],color='#246cab',lw=.65,label='u₀')
    line2,=ax.plot([],[],[],color='#dc6848',lw=.65,label='u₀+δ');ax.legend(loc='upper right',fontsize=8)
    point1,=ax.plot([],[],[],'o',color='#246cab',ms=3);point2,=ax.plot([],[],[],'o',color='#dc6848',ms=3)
    text=ax.text2D(.02,.92,'',transform=ax.transAxes,fontsize=9)
    frames=np.unique(np.linspace(1,len(t),c['animation_frames']).astype(int))
    def update(end):
        for line,point,values in [(line1,point1,y),(line2,point2,yp)]:
            line.set_data_3d(*values[:end].T);point.set_data_3d(*values[end-1:end].T)
        text.set_text(f"t={t[end-1]:.2f}; d={data['distance'][end-1]:.3e}")
        return line1,line2,point1,point2,text
    animation=FuncAnimation(fig,update,frames=frames,blit=False)
    animation.save(out/'animation.gif',writer=PillowWriter(fps=c['animation_fps']),dpi=90)
    plt.close(fig)
