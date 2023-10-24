#This code was not written by me. It was writen by Arefe Ghazi Nezami

import numpy as np
import dedalus.public as d3
from dedalus.tools.logging import *
logger = logging.getLogger(__name__) #what is this?


# Parameters
Lx, Ly, Lz = 2*np.pi, 2*np.pi, 2*np.pi
Nx, Ny, Nz = 256, 256, 256
dealias = 3/2
N02 = 1
om = 0.5
zf = np.pi
zd = 0.75
sigma = 0.25
wsigma = 0.25
T_d = np.pi/N02
zw = 1.5
P = 2*np.pi/om
pw = 1
nu = kappa = 5e-4
c = 0.1
Pw = 1.5
vol = Lx*Ly*Lz
amp = 0.0025
kx = 8*np.pi/Lx
Tw = 3*P



timestepper = d3.SBDF2
max_timestep = 2*np.pi/np.sqrt(N02)*0.25
stop_sim_time = 801
dtype = np.float64

#Bases
coords = d3.CartesianCoordinates('x', 'y', 'z')
dist = d3.Distributor(coords, dtype=dtype , mesh=[32,16])
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
x, y, z = dist.local_grids(xbasis, ybasis, zbasis)

#Fields
p = dist.Field(name='p', bases=(xbasis,ybasis,zbasis))
b = dist.Field(name='b', bases=(xbasis,ybasis,zbasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis,zbasis))
tau_p = dist.Field(name='tau_p')
tau_b1 = dist.Field(name='tau_b1', bases=(xbasis , ybasis))
tau_b2 = dist.Field(name='tau_b2', bases=(xbasis, ybasis))
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=(xbasis,ybasis))
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=(xbasis , ybasis))
t = dist.Field(name='t')


sinx = dist.Field(name='sin(x)', bases=xbasis)
sinx['g'] = np.sin(kx*x)
sinx = d3.Grid(sinx)

cosx = dist.Field(name='cos(x)', bases=xbasis)
cosx['g'] = np.cos(kx*x)
cosx = d3.Grid(cosx)


#damping
damping = dist.Field(name='D', bases=zbasis)
damping['g'] = (np.exp(-((3.5* (z-zd)/zd)**2/2)))

# Wave forcing
wforcing = dist.Field(name='F', bases=zbasis)
wforcing['g'] = amp*np.exp( -(z-zw)**2/2/wsigma**2 )
wforcing = d3.Grid(wforcing)



#Forcing parameters
kx = xbasis.local_elements()[0]
ky = ybasis.local_elements()[0]
f2D = dist.VectorField(coords, name= 'f2D' , bases=(xbasis, ybasis, zbasis))
f3D = dist.VectorField(coords, name= 'f3D' , bases=(xbasis, ybasis, zbasis))
f2D.preset_scales(dealias)
f3D.preset_scales(dealias)
ky, kx= np.meshgrid(ky,kx)
kh = np.sqrt(kx**2+ky**2)

forcing = dist.Field(name='F', bases=(xbasis, ybasis, zbasis))
forcing['g'] = 0.5*(1+ np.tanh(2*(z-zf)))
forcing.change_scales(dealias)
forcing_shape = forcing['g'].shape
forcing_grid = forcing['g'].reshape((1, *forcing_shape))
#forcing = d3.Grid(forcing)


mask = (kh > 2.5) & (kh < 3.5)
knum = np.sum(mask)

rand = np.random.RandomState(seed=34713)


#
ex, ey, ez = coords.unit_vector_fields(dist)
lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)
grad_u = d3.grad(u) + ez*lift(tau_u1) # First-order reductio
grad_b = d3.grad(b) + ez*lift(tau_b1) # First-order reductio
stress = 1/2*(d3.grad(u) + d3.grad(u).T)






# Problem

problem = d3.IVP([p, b, u, tau_p, tau_b1, tau_b2, tau_u1, tau_u2], time=t, namespace=locals())
problem.add_equation("trace(grad_u) + tau_p = 0")
problem.add_equation("dt(b) - kappa*div(grad_b) + lift(tau_b2) +N02*(ez@u) = - u@grad(b) - damping*b + wforcing*(np.cos(om*t)*cosx+ np.sin(om*t)*sinx)*0.5*(1+np.tanh( 2*(t-Tw)/Tw ))")
problem.add_equation("dt(u) - nu*div(grad_u) + grad(p) - b*ez + lift(tau_u2) = "  +
" f3D - u@grad(u) - damping*u -ex*wforcing*np.sqrt(N02**2/om**2-1)*(om/N02**2)*(np.cos(om*t)*sinx - np.sin(om*t)*cosx)*0.5*(1+np.tanh( 2*(t-Tw)/Tw ))" +
                    "-ez*(wforcing*(om/N02**2)*(np.cos(om*t)*sinx - np.sin(om*t)*cosx))*0.5*(1+np.tanh(2*(t-Tw)/Tw ))" )
problem.add_equation("b(z=0) = 0")
problem.add_equation("ez@u(z=0) = 0")
problem.add_equation("ex@stress(z=0)@ez=0")
problem.add_equation("ey@stress(z=0)@ez=0")
#problem.add_equation("ex@(ez@stress(z=0))=0")
problem.add_equation("b(z=Lz) = 0")
problem.add_equation("ez@u(z=Lz) = 0")
problem.add_equation("ex@stress(z=Lz)@ez=0")
problem.add_equation("ey@stress(z=Lz)@ez=0")
problem.add_equation("integ(p) = 0") # Pressure gauge



# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

epsilon = nu*d3.curl(u)@d3.curl(u)
uh = (u)@(u)

# IC
#u.fill_random('g', seed=4852, distribution='normal', scale=5e-4)
#u['g'] *= 0.5*(1+ np.tanh(10*(z-zf)))

initial_timestep = solver.load_state('checkpoints_s7.h5')
file_handler_mode = 'append'
#Analysis
slices = solver.evaluator.add_file_handler('output/slices', sim_dt=1, max_writes=50)
slices.add_task(b(x=np.pi), name='buoyancy_x_mid')
slices.add_task(epsilon(x=np.pi), name='epsilon_x_mid')
slices.add_task(u(x=np.pi), name= 'U_x_mid')
slices.add_task(p(x=np.pi), name= 'p_x_mid')
slices.add_task(b(y=np.pi), name='buoyancy_y_mid')
slices.add_task(epsilon(y=np.pi), name='epsilon_y_mid')
slices.add_task(u(y=np.pi), name= 'U_y_mid')
slices.add_task(p(y=np.pi), name= 'p_y_mid')
slices.add_task(b(z=1.5*np.pi), name='buoyancy_z_midt')
slices.add_task(epsilon(z=1.5*np.pi), name='epsilon_z_midt')
slices.add_task(u(z=1.5*np.pi), name= 'U_z_midt')
slices.add_task(b(z=0.75*np.pi), name='buoyancy_z_midnt')
slices.add_task(epsilon(z=0.75*np.pi), name='epsilon_z_midnt')
slices.add_task(u(z=0.75*np.pi), name= 'U_z_midnt')
slices.add_task(b(z=np.pi), name='buoyancy_z_t')
slices.add_task(epsilon(z=np.pi), name='epsilon_z_t')
slices.add_task(u(z=np.pi), name= 'U_z_t')
slices.add_task(d3.ave(b,'y'), name= 'b_yave')


checkpoints = solver.evaluator.add_file_handler('output/checkpoints' , wall_dt=59*60 , max_writes=1)
checkpoints.add_tasks(solver.state)

scalars = solver.evaluator.add_file_handler('output/scalars' , sim_dt=0.1)
scalars.add_task(d3.integ(1/2*u@u) , name='KE')
scalars.add_task(d3.integ(epsilon)/vol , name='dissipation')
scalars.add_task(d3.integ(u@f3D) , name='power')
scalars.add_task(d3.ave(ez@u) , name='Wbar')
scalars.add_task(np.sqrt(d3.ave(ex@u**2) + d3.ave(ey@u**2)) , name= 'uh')
scalars.add_task(np.sqrt(d3.integ(ex@u**2)) , name= 'urms')
scalars.add_task(d3.integ(ey@u**2) , name= 'vrms')
scalars.add_task(d3.integ(ez@u**2) , name= 'wrms')
scalars.add_task(d3.integ(ex@u*ey@u) , name= 'uv')
scalars.add_task(d3.integ(ex@u*ez@u) , name= 'uw')
scalars.add_task(d3.integ(ey@u*ez@u) , name= 'vw')
scalars.add_task(d3.integ(ez@u*b) , name='bflux')



# CFL
CFL = d3.CFL(solver, initial_dt=0.01*max_timestep, cadence=1, safety=0.2, threshold=0.05,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=1)
flow.add_property(np.sqrt(u@u)/nu, name='Re')

Puf_op = d3.integ(u@f3D)
Pff_op = d3.integ(f3D@f3D)

# Main loop
startup_iter = 10
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        phases = rand.uniform(0, 2 * np.pi, knum)
        f2D.preset_layout('c')
        f2D['c'] *= 0
        f2D['c'][0, :, :, 0][mask] = (ky / kh)[mask] * np.cos(phases) / (np.sqrt(np.pi * kh[mask]))/np.sqrt(timestep)
        f2D['c'][1, :, :, 0][mask] = (-kx / kh)[mask] * np.cos(phases) / (np.sqrt(np.pi * kh[mask]))/np.sqrt(timestep)
        f3D.preset_scales(dealias)
        f3D.preset_layout('g')
        f3D['g'] = f2D['g']*forcing_grid
        Puf = flow.reducer.global_max(Puf_op.evaluate()['g'])
        Pff = 0.5*timestep*flow.reducer.global_max(Pff_op.evaluate()['g'])
        c = min(np.abs((-Puf+np.sqrt(Puf**2 + 4*Pff*Pw))/(2*Pff)) , np.abs((-Puf-np.sqrt(Puf**2 + 4*Pff*Pw))/(2*Pff)))
        f3D['c'] *= c
        #f['g'] *= forcing
        solver.step(timestep)
        if (solver.iteration-1) % 1 == 0:
            max_Re = flow.max('Re')
            logger.info('Iteration=%i, Time=%e, dt=%e, max(Re)=%f' %(solver.iteration, solver.sim_time, timestep, max_Re))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
#finally:
#    solver.log_stats()

