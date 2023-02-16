import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p


from design_opt import *
import aerosandbox.numpy as np

try:
    sol = opti.solve()
except RuntimeError:
    sol = opti.debug
s = lambda x: sol.value(x)

airplane.substitute_solution(sol)
cruise_op_point.substitute_solution(sol)
# dyn.substitute_solution(sol)
mass_props_TOGW.substitute_solution(sol)

opti2 = asb.Opti()

vels = np.linspace(5, 25, 400)

op_point_sweep = asb.OperatingPoint(
    atmosphere=asb.Atmosphere(altitude=0),
    velocity = vels,
    alpha=opti2.variable(
        init_guess=5 * np.ones_like(vels),
        lower_bound=-10,
        upper_bound=15
    )
)
aero_sweep = asb.AeroBuildup(
    airplane=airplane,
    op_point=op_point_sweep,
    xyz_ref=mass_props_TOGW.xyz_cg
).run()

aero_sweep['D'] = aero_sweep['D'] + s(
        fuse_cfd_results['CDA'] * cruise_op_point.dynamic_pressure() +  # Fuselage
        (  # Sponson
                2 *
                lib_aero.Cf_flat_plate(cruise_op_point.reynolds(reference_length=wing_sponson_length)) *
                1.47 *  # Form factor, calibrated to fuselage shape
                wing_sponson_wetted_area *
                cruise_op_point.dynamic_pressure() *
                3  # rough WAG to account for sponson mounts
        ) +
        (  # Boom
                lib_aero.Cf_flat_plate(cruise_op_point.reynolds(reference_length=boom_length)) *
                boom_wetted_area *
                cruise_op_point.dynamic_pressure()
        )
)

lift_residual = aero_sweep['L'] - mass_props_TOGW.mass * 9.81

opti2.minimize(np.sum(lift_residual ** 2))

sol2 = opti2.solve()

LDs = sol2.value(aero_sweep['L'] / aero_sweep['D'])

trimmable_points = (np.abs(sol2.value(lift_residual)) < 1e-3)

sink_rates = sol2.value(vels / LDs)

fig, ax = plt.subplots(
    figsize=(5,4.5)
)
line = plt.plot(
    vels[trimmable_points],
    -sink_rates[trimmable_points],
    label="Ideal Piloting"
)
plt.plot(
    vels[trimmable_points],
    -sink_rates[trimmable_points] / 0.75,
    # "--",
    # color=line[0].get_color(),
    label="Realistic Piloting\n(RC, with wind)"
)
plt.plot(
    [9.38],
    [-0.63],
    ".k",
)
plt.annotate(
            text="SEAWAY-Mini\nCruise",
            xy=(9.38, -0.63),
            xytext=(7, -1.5),
            xycoords="data",
            ha='center',
            arrowprops={
                "color"     : "k",
                "width"     : 0.25,
                "headwidth" : 4,
                "headlength": 6,
                "shrink"    : 0.1
            }
        )

plt.xlim(left=0)
plt.ylim(top=0)
p.set_ticks(5, 1, 0.5, 0.1)
p.show_plot(
    "Glide Polar for SEAWAY-Mini",
    "Airspeed [m/s]",
    "Sink Rate [m/s]"
)