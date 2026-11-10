"""Проверки опираются на точные решения, инварианты и независимый интегратор."""
import unittest
import numpy as np
from scipy.integrate import solve_ivp
from lorenz import rhs,jacobian,equilibria,Parameters
from numerical_methods import euler_step,rk4_step,integrate,time_grid

class NumericsTests(unittest.TestCase):
    def test_rhs_known_value(self):
        np.testing.assert_allclose(rhs(7,np.ones(3)),[0,26,-5/3],rtol=0,atol=1e-14)

    def test_euler_known_step(self):
        np.testing.assert_allclose(euler_step(rhs,0,np.ones(3),.01),[1,1.26,1-1/60],rtol=0,atol=1e-14)

    def test_rk4_exponential(self):
        sol=integrate(lambda t,u:u,[1],0,1,.01)
        self.assertLess(abs(sol.y[-1,0]-np.e),3e-10)

    def test_rk4_nonautonomous(self):
        # u'=t*u, u(0)=1 => u(t)=exp(t²/2), проверяет времена стадий.
        sol=integrate(lambda t,u:t*u,[1],0,1,.01)
        self.assertLess(abs(sol.y[-1,0]-np.exp(.5)),1e-8)

    def test_rk4_order(self):
        errors=[abs(integrate(lambda t,u:u,[1],0,1,h).y[-1,0]-np.e) for h in [.1,.05,.025]]
        for a,b in zip(errors,errors[1:]):self.assertTrue(3.8<np.log2(a/b)<4.2)

    def test_euler_order(self):
        errors=[abs(integrate(lambda t,u:u,[1],0,1,h,'euler').y[-1,0]-np.e) for h in [.02,.01,.005]]
        for a,b in zip(errors,errors[1:]):self.assertTrue(.9<np.log2(a/b)<1.1)

    def test_short_agreement(self):
        rk=integrate(rhs,[1,1,1],0,2,.001)
        ref=solve_ivp(rhs,(0,2),[1,1,1],method='DOP853',rtol=3e-14,atol=1e-15,t_eval=rk.t)
        self.assertTrue(ref.success)
        self.assertLess(np.linalg.norm(rk.y-ref.y.T,axis=1).max(),2e-6)

    def test_equilibria_and_spectra(self):
        eq=equilibria();self.assertEqual(len(eq),3)
        for u in eq.values():np.testing.assert_allclose(rhs(0,u),0,rtol=0,atol=5e-14)
        self.assertTrue(np.any(np.real(np.linalg.eigvals(jacobian(eq['C+'])))>0))
        self.assertEqual(len(equilibria(Parameters(rho=1))),1)
        self.assertEqual(len(equilibria(Parameters(rho=.5))),1)

    def test_jacobian_independent_difference(self):
        u=np.array([2.3,-1.4,8.1]);eps=1e-5
        fd=np.column_stack([(rhs(0,u+eps*e)-rhs(0,u-eps*e))/(2*eps) for e in np.eye(3)])
        np.testing.assert_allclose(jacobian(u),fd,rtol=0,atol=1e-8)

    def test_grid_and_evaluation_counts(self):
        t=time_grid(0,40,.01);self.assertEqual(len(t),4001);self.assertEqual(t[-1],40)
        for method,k in [('rk4',4),('euler',1)]:
            calls=[]
            def fun(t,u):calls.append(t);return -u
            sol=integrate(fun,[1],0,1,.1,method)
            self.assertEqual(sol.nfev,len(calls));self.assertEqual(len(calls),10*k)
        for h in [0,-1,.03]:
            with self.assertRaises(ValueError):time_grid(0,1,h)

    def test_standard_trajectories_finite(self):
        for method in ['euler','rk4']:
            for u in [[1,1,1],[1+1e-8,1,1]]:
                sol=integrate(rhs,u,0,40,.01,method)
                self.assertTrue(np.isfinite(sol.y).all())

    def test_overflow_not_silent(self):
        with self.assertRaises(FloatingPointError):
            integrate(lambda t,u:u*u,[1e200],0,1,.1)

if __name__=='__main__':unittest.main()
