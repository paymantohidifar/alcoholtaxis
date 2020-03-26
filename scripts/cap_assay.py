from fenics import *
import numpy as np
import matplotlib.pyplot as plt


def integrator(u, dx):
    """Take solution vector (u) and domain (dx) and
    return the integration value"""
    return assemble(u*dx)


def make_figures(z, t, u_states, save_path):

    try:
        # concentration profiles along center line at time = 0, 10, and 30 min
        plt.figure(1)
        plt.plot(z, u_states[0].get('u_z'), 'k-')
        plt.plot(z, u_states[599].get('u_z'), 'b-')
        plt.plot(z, u_states[1799].get('u_z'), 'r-')
        plt.title('Ethanol concentration profile along capillary at different time points')
        plt.xlabel('Position (mm)')
        plt.ylabel('Ethanol concentration (mM)')
        plt.legend(['t = 0', 't = 10 min', 't = 30 min'])
        plt.savefig(save_path + 'fig_1.eps', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                    format='eps', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    except IndexError:
        print('Cannot generate figure 1!')

    # Concentration profiles over 30 min at z = 0.1, 0.25, and 0.5 mm along center line on z axis
    plt.figure(2)
    plt.semilogx(t, [state['u_z'][state['z_pos'].index(0.1)] for state in u_states], 'k-')
    plt.semilogx(t, [state['u_z'][state['z_pos'].index(0.25)] for state in u_states], 'b-')
    plt.semilogx(t, [state['u_z'][state['z_pos'].index(0.5)] for state in u_states], 'r-')
    plt.title('Ethanol concentration evolution in time at different points along z axis')
    plt.xlabel('Time (min)')
    plt.ylabel('Ethanol concentration (mM)')
    plt.legend(['z = 0.1 mm', 'z = 0.25 mm', 'z = 0.5 mm'])
    plt.savefig(save_path + 'fig_2.eps', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                format='eps', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

    # Total mass change in capillary, pond, and whole system (control simulation for mass conservation)
    plt.figure(3)
    plt.plot(t, [state['cap_mass'] for state in u_states], 'k-')
    plt.plot(t, [state['pond_mass'] for state in u_states], 'b-')
    plt.plot(t, [state['cap_mass'] + state['pond_mass'] for state in u_states], 'r-')
    plt.title('Total mass change in system')
    plt.xlabel('Time (min)')
    plt.ylabel('Ethanol mass')
    plt.legend(['Capillary', 'Pond', 'Total'])

    plt.savefig(save_path + 'fig_3.eps', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                format='eps', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    plt.show()


def write_data(z, t, u_states, write_path):
    try:
        with open(write_path + 'fig1.txt', 'w') as file:
            file.write('\t'.join(['Position', 't=0', 't=10', 't=30']) + '\n')
            file.writelines('\t'.join([str(z[n]), str(u_states[0].get('u_z')[n]), str(u_states[599].get('u_z')[n]),
                                       str(u_states[1799].get('u_z')[n])]) + '\n' for n in range(0, len(z)))
    except IndexError:
        print('Cannot write data for figure 1!')

    with open(write_path + 'fig2.txt', 'w') as file:
        file.write('\t'.join(['Time', 'z=0.1', 'z=0.25', 'z=0.5']) + '\n')
        file.writelines('\t'.join([str(t[n]), str(state['u_z'][state['z_pos'].index(0.1)]),
                                   str(state['u_z'][state['z_pos'].index(0.25)]),
                                   str(state['u_z'][state['z_pos'].index(0.5)])]) + '\n'
                        for n, state in enumerate(u_states))


def main():
    """ solves spatio-temporal dynamics of ethanol concentration
    in capillary assay using finite-element method """

    # Only verbose errors
    set_log_level(40,)

    # set finite element positional tolerance
    tol = 1e-14

    # Create domain mesh
    mesh = Mesh("../mesh/cap3d.xml")

    # Define capillary and pond sub-domains
    subdomain_cap = CompiledSubDomain('x[2] <= 0. + tol', tol=tol)
    subdomain_pond = CompiledSubDomain('x[2] >= 0. - tol', tol=tol)
    mf = MeshFunction('size_t', mesh, mesh.topology().dim())
    subdomain_cap.mark(mf, 0)
    subdomain_pond.mark(mf, 1)
    dx = Measure('dx', domain=mesh, subdomain_data=mf)

    # Define 'Lagrange' function space using Continuous Galerkin (CG) elements
    V = FunctionSpace(mesh, 'CG', 1)

    # Set up problem parameters
    u_ic = Expression('C0*(x[2] < 0.)', degree=0, C0=50.)
    u0 = interpolate(u_ic, V)

    # Define natural Neumann boundary condition (g = 0) with no source (f = 0)
    g = Constant(0)
    f = Constant(0)

    # Define diffusion coefficient of ethanol in water (mm2/s)
    D = Constant(1.23e-3)
    # print('Ethanol diffusion coefficient: {} mm2/s'.format(D.values()[0]))

    # Define time variables
    ftime = 1800         # 30-min long
    num_steps = 1800
    dt = ftime / num_steps  # dt = 1 s

    # Define variational problem (Euler backward implicit method is used for time integration
    # (L-stable: independent of dt size))
    u = TrialFunction(V)
    v = TestFunction(V)
    a = v*u*dx + dt*D*dot(grad(v), grad(u))*dx
    L = v*u0*dx + dt*D*v*g*ds + dt*v*f*dx

    # Time-stepping
    u = Function(V)
    t = 0
    time_series = []
    u_states = []
    z_axis = np.arange(-5, 2, 0.01)
    for n in range(num_steps+1):

        # perform evaluations every 1 s
        if np.mod(n, 1) == 0:
            time_series.append(t)

            # Compute concentration along center line z axis
            u_state = {'z_pos': [round(z, 2) for z in z_axis],
                       'u_z': [u0(0, 0, z) for z in z_axis],
                       'cap_mass': integrator(u0, dx(0)),
                       'pond_mass': integrator(u0, dx(1))}
            u_states.append(u_state)

            # Verbose progress
            print('Progress: {:.2f}%'.format(t / ftime * 100))

        # Compute new solutions and update
        t += dt
        solve(a == L, u)
        u0.assign(u)

    # Output solutions in text files
    resultsdir = '../results/'
    write_data(z_axis, time_series, u_states, resultsdir)

    # generate figures
    figuredir = '../results/figures/'
    make_figures(z_axis, time_series, u_states, figuredir)


# ------------------------------
if __name__ == '__main__':
    print('PDE solver begins:\n' + 100 * '-')
    main()
    print('End\n' + 100 * '-')
