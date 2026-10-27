"""Independent physical checks of open-system, moving-CG dynamics."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from VehicleData import build_falcon9
from dynamics import State, Loads, rk4_step, accelerations, acceleration_truth
from mass_properties import MassProperties
from mass_flow import ExhaustPort, uniform_disk_gyration
from propagation import propagate
from frames import body_to_initial_body_matrix


class LumpedRocket:
    mass_model_name = 'analytic_fixture'
    def propellant_bounds(self,configuration):
        return np.zeros(2),np.array([100.,0.])
    def mass_properties(self,configuration,propellant_kg,consumption_kgps=(0,0),consumption_rate_kgps2=(0,0)):
        m=100+propellant_kg[0]
        return MassProperties(m,np.zeros(3),np.eye(3)*m,np.eye(3)*(-consumption_kgps[0]),
                               -consumption_kgps[0],np.zeros(3),(),-consumption_rate_kgps2[0],np.zeros(3))


def centered_provider(flow,flow_change=0,thrust=0):
    return lambda t,s: Loads([thrust,0,0],[0,0,0],[0,0,0],[flow,0],
                              (ExhaustPort(0,[0,0,0],flow,flow_change),))


class PoweredDynamicsChecks(unittest.TestCase):
    def test_tsiolkovsky_velocity_and_position(self):
        state=State([0,0,0],[0,0,0],[1,0,0,0],[0,0,0],[100,0])
        q,ve,duration=2.,1000.,10.
        result=propagate(0,state,duration,.05,LumpedRocket(),None,centered_provider(q,thrust=q*ve))
        final=result.states[-1]; m0,m1=200.,180.
        self.assertAlmostEqual(final.velocity_body_mps[0],ve*np.log(m0/m1),places=9)
        self.assertAlmostEqual(final.position_body_m[0],ve*(duration-m1/q*np.log(m0/m1)),places=8)

    def test_missing_or_unbalanced_ports_are_rejected(self):
        with self.assertRaises(ValueError): Loads([0,0,0],[0,0,0],[0,0,0],[1,0])
        with self.assertRaisesRegex(ValueError,'balanced'):
            Loads([0,0,0],[0,0,0],[0,0,0],[1,0],(ExhaustPort(0,[0,1,0],1),))

    def test_analytic_cg_derivatives_match_smooth_inventory(self):
        model=build_falcon9().mass_model()
        prop=np.array([200000.,107500.]); q=np.array([1500.,0]); dq=np.array([30.,0]); h=.002
        now=model.mass_properties(None,prop,q,consumption_rate_kgps2=dq)
        plus=model.mass_properties(None,prop-q*h-.5*dq*h*h)
        minus=model.mass_properties(None,prop+q*h-.5*dq*h*h)
        np.testing.assert_allclose((plus.cg_station_body_m-minus.cg_station_body_m)/(2*h),now.cg_rate_body_mps,rtol=1e-8,atol=1e-10)
        np.testing.assert_allclose((plus.cg_station_body_m-2*now.cg_station_body_m+minus.cg_station_body_m)/(h*h),
                                   now.cg_acceleration_body_mps2,rtol=2e-4,atol=2e-8)

    def test_variable_flow_has_nonzero_second_derivative(self):
        model=build_falcon9().mass_model()
        p=model.mass_properties(None,[300000,107500],[1000,0],consumption_rate_kgps2=[100,0])
        self.assertAlmostEqual(p.mass_acceleration_kgps2,-100)
        self.assertNotEqual(p.cg_acceleration_body_mps2[0],0)

    def test_roll_angular_momentum_leaves_with_the_jet(self):
        model=LumpedRocket()
        state=State([0,0,0],[0,0,0],[1,0,0,0],[.2,0,0],[100,0])
        # Exit axial gyration=1 equals fixture Jxx/m=1; spin rate stays constant.
        loads=Loads([0,0,0],[0,0,0],[0,0,0],[2,0],
                    (ExhaustPort(0,[0,0,0],2,gyration_body_m2=uniform_disk_gyration(np.sqrt(2))),))
        _,alpha=accelerations(state,loads,model.mass_properties(None,state.propellant_kg,[2,0]))
        np.testing.assert_allclose(alpha,[0,0,0],atol=1e-14)

    def test_cg_equations_match_fixed_reference_coupled_balance(self):
        v=build_falcon9(); model=v.mass_model(); state=v.to_state(); omega=np.array([.01,.02,-.015])
        state=State(state.position_body_m,state.velocity_body_mps,state.quaternion_body_to_initial_body_wxyz,omega,state.propellant_kg)
        loads=centered_provider(2000,40,6e6)(0,state)
        p=model.mass_properties(None,state.propellant_kg,[2000,0],consumption_rate_kgps2=[40,0])
        a,alpha=accelerations(state,loads,p)
        c,cp,cpp=p.cg_station_body_m,p.cg_rate_body_mps,p.cg_acceleration_body_mps2
        m,md,mdd=p.mass_kg,p.mass_rate_kgps,p.mass_acceleration_kgps2
        cross=lambda r:np.array([[0,-r[2],r[1]],[r[2],0,-r[0]],[-r[1],r[0],0]])
        S=lambda r:np.dot(r,r)*np.eye(3)-np.outer(r,r)
        JO=p.inertia_body_kgm2+m*S(c)
        JOdot=p.inertia_rate_body_kgm2ps+md*S(c)+m*(2*np.dot(c,cp)*np.eye(3)-np.outer(cp,c)-np.outer(c,cp))
        F=loads.force_body_n
        linear=F-m*np.cross(omega,np.cross(omega,c))-2*m*np.cross(omega,cp)-m*cpp-2*md*cp+2*md*np.cross(omega,-c)-mdd*c
        angular=loads.moment_body_nm+np.cross(c,F)-np.cross(omega,JO@omega)-JOdot@omega
        block=np.block([[m*np.eye(3),-m*cross(c)],[m*cross(c),JO]])
        result=np.linalg.solve(block,np.r_[linear,angular])
        cg_a=result[:3]+np.cross(result[3:],c)+np.cross(omega,np.cross(omega,c))+2*np.cross(omega,cp)+cpp
        np.testing.assert_allclose(alpha,result[3:],rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(a,cg_a,rtol=1e-10,atol=1e-12)

    def test_state_round_trip_with_cg_migration(self):
        v=build_falcon9(); v.omega=np.array([.02,-.01,.03])
        v.position=np.array([10.,20.,30.]); v.velocity=np.array([2.,3.,4.])
        x,u=v.position.copy(),v.velocity.copy()
        state=v.to_state(consumption_by_stage={'stage1':1000},consumption_derivative_by_stage={'stage1':20})
        v.apply_state(state,5,consumption_by_stage={'stage1':1000},consumption_derivative_by_stage={'stage1':20})
        np.testing.assert_allclose(v.position,x,atol=1e-12)
        np.testing.assert_allclose(v.velocity,u,atol=1e-12)

    def test_fuel_event_stops_at_limiting_inventory(self):
        lower,_=build_falcon9().separate()
        minimum=lower.mass_model().propellant_bounds(None)[0][0]
        lower.stages[0].propellant_mass=minimum+10
        state=lower.to_state(consumption_by_stage={'stage1':2})
        result=propagate(0,state,10,1,lower.mass_model(),None,centered_provider(2,thrust=6000))
        self.assertEqual(result.stop_reason,'propellant_limit')
        self.assertAlmostEqual(result.times_s[-1],5,places=6)
        self.assertGreaterEqual(result.states[-1].propellant_kg[0],minimum-1e-7)

    def test_sensor_acceleration_matches_base_velocity_change(self):
        v=build_falcon9(); model=v.mass_model()
        state=v.to_state(consumption_by_stage={'stage1':1000})
        loads=centered_provider(1000,thrust=3e6)(0,state)
        p=model.mass_properties(None,state.propellant_kg,[1000,0])
        truth=acceleration_truth(state,loads,p,sensor_station_body_m=[0,0,0])
        h=1e-3; plus=rk4_step(0,state,h,model,None,centered_provider(1000,thrust=3e6))
        pplus=model.mass_properties(None,plus.propellant_kg,[1000,0])
        R=body_to_initial_body_matrix(state.quaternion_body_to_initial_body_wxyz)
        Rp=body_to_initial_body_matrix(plus.quaternion_body_to_initial_body_wxyz)
        u=R@(state.velocity_body_mps-p.cg_rate_body_mps)
        up=Rp@(plus.velocity_body_mps-pplus.cg_rate_body_mps)
        np.testing.assert_allclose((up-u)/h,truth['sensor_acceleration_body_mps2'],rtol=1e-5,atol=2e-5)

    def test_rising_flow_fuel_event_is_resolved_without_negative_inventory(self):
        lower,_=build_falcon9().separate()
        minimum=lower.mass_model().propellant_bounds(None)[0][0]
        lower.stages[0].propellant_mass=minimum+10
        state=lower.to_state(consumption_by_stage={'stage1':1},consumption_derivative_by_stage={'stage1':1})
        provider=lambda t,s: centered_provider(1+t,1,3000*(1+t))(t,s)
        result=propagate(0,state,10,1,lower.mass_model(),None,provider)
        self.assertEqual(result.stop_reason,'propellant_limit')
        self.assertAlmostEqual(result.times_s[-1],np.sqrt(21)-1,places=6)
        remaining=result.states[-1].propellant_kg[0]-minimum
        self.assertGreaterEqual(remaining,-1e-7)
        self.assertLessEqual(remaining,1e-6)


if __name__=='__main__': unittest.main()
