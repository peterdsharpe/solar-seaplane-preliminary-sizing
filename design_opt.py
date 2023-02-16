import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.library.aerodynamics as lib_aero
from aerosandbox.library import mass_structural as lib_mass_struct
from aerosandbox.library import power_solar as lib_solar
from aerosandbox.library import propulsion_electric as lib_prop_elec
from aerosandbox.library import propulsion_propeller as lib_prop_prop
import aerosandbox.tools.units as u
import copy
from pathlib import Path

opti = asb.Opti(
    freeze_style='float',
    # variable_categories_to_freeze='all'
)

make_plots = True

##### Section: Parameters

n_panels_spanwise_inboard = opti.parameter(6)  # number of panels per side
n_panels_chordwise_inboard = opti.parameter(3)
n_panels_spanwise_outboard = opti.parameter(9)  # number of panels per side
n_panels_chordwise_outboard = opti.parameter(2)

panel_spacing = 0.127 - 1e-6  # center to center, the closest you can get
solar_area_per_panel = (0.125) ** 2

wing_extra_span = 1.5 * u.inch

wing_span = (
        2 * (n_panels_spanwise_inboard + n_panels_spanwise_outboard) * panel_spacing +
        wing_extra_span
)

panel_chord_coverage_fraction = 0.90

wing_y_break_fraction = n_panels_spanwise_inboard / (n_panels_spanwise_inboard + n_panels_spanwise_outboard)

battery_specific_energy_Wh_kg = 150  # at pack level
battery_voltage = 4 * 3.7
battery_reserve_time = 20 * u.minute

solar_cell_efficiency = 0.243 * 0.9  # Sunpower
rho_solar_cells = 0.425 * 1.1 * 1.15  # kg/m^2, solar cell area density. Sunpower.
latitude = 42.36  # Boston
day_of_year = 91  # April 1, roughly
time_after_solar_noon = 2 * u.hour

wing_dihedral_angle_deg = 6

##### Section: Vehicle Overall Specs

cruise_op_point = asb.OperatingPoint(
    velocity=opti.variable(
        init_guess=14,
        lower_bound=4,
        upper_bound=50,
        log_transform=True
    ),
    alpha=opti.variable(
        init_guess=0,
        lower_bound=-10,
        upper_bound=10,
    )
)

design_mass_TOGW = opti.variable(
    init_guess=10,
    lower_bound=1e-3,
)
design_mass_TOGW = np.maximum(design_mass_TOGW, 1e-3)

LD_cruise = opti.variable(
    init_guess=15,
    lower_bound=0.1,
    log_transform=True,
)  # TODO close this loop

g = 9.81

target_climb_angle = 25  # degrees

design_thrust_cruise_total = (
        design_mass_TOGW * g / LD_cruise  # cruise component
)
design_thrust_climb_total = (
        design_mass_TOGW * g / LD_cruise +
        design_mass_TOGW * g * np.sind(target_climb_angle)
)

# TODO compute takeoff thrust needs


##### Section: Vehicle Definition

"""
Coordinate system:

Geometry axes. Datum (0, 0, 0) is coincident with the quarter-chord-point of the centerline cross section of the main 
wing.
    
"""

### Basic configuration
x_tail = 0 + wing_span * 0.38

### Wing
wing_root_chord = (n_panels_chordwise_inboard * panel_spacing) / 0.90
wing_tip_chord = (n_panels_chordwise_outboard * panel_spacing + 1 * u.inch) / 0.90

airfoils = {
    name: asb.Airfoil(
        name=name,
        coordinates=Path(".") / "airfoils" / f"{name}.dat"
    )
    for name in [
        "ag34",
        "ag35",
        "ag36",
        "ht14",
    ]
}

for v in airfoils.values():
    v.generate_polars(
        cache_filename=f"cache/{v.name}.json",
        alphas=np.linspace(-12, 12, 30)
    )


def wing_rot(xyz):
    dihedral_rot = np.rotation_matrix_3D(
        angle=np.radians(wing_dihedral_angle_deg),
        axis="X"
    )

    return dihedral_rot @ np.array(xyz)


wing = asb.Wing(
    name="Wing",
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=wing_rot([
                -wing_root_chord / 4,
                0,
                0,
            ]),
            chord=wing_root_chord,
            twist=4.70,
            airfoil=airfoils['ag34']
        ),
        asb.WingXSec(
            xyz_le=wing_rot([
                -wing_root_chord / 4,
                wing_span / 2 * wing_y_break_fraction,
                0
            ]),
            chord=wing_root_chord,
            twist=3.50,
            airfoil=airfoils['ag34']
        ),
        asb.WingXSec(
            xyz_le=wing_rot([
                -wing_root_chord / 4,
                wing_span / 2,
                0
            ]),
            chord=wing_tip_chord,
            twist=3.00,
            airfoil=airfoils['ag36']
        )
    ]
)

# Compute aileron properties
aileron_inboard_chord = wing_root_chord * 0.90 - n_panels_chordwise_outboard * panel_spacing
aileron_outboard_chord = wing_tip_chord * 0.90 - n_panels_chordwise_outboard * panel_spacing
aileron_total_area = 2 * (
        (aileron_inboard_chord + aileron_outboard_chord) / 2 *
        (1 - wing_y_break_fraction) * wing_span / 2
)
aileron_area_fraction = aileron_total_area / wing.area()

### VStab
vstab_volume_coefficient = 0.030  # Potentially as low as 0.0111
vstab_aspect_ratio = 2.5
vstab_taper_ratio = 0.7

vstab_area = vstab_volume_coefficient * wing.area() * wing.span() / x_tail
vstab_span = (vstab_area * vstab_aspect_ratio) ** 0.5
vstab_chord = (vstab_area / vstab_aspect_ratio) ** 0.5
vstab_airfoil = airfoils['ht14']

vstab = asb.Wing(
    name="VStab",
    xsecs=[
        asb.WingXSec(
            xyz_le=[
                -vstab_chord * vstab_taper_ratio ** -0.5,
                0,
                0
            ],
            chord=vstab_chord * vstab_taper_ratio ** -0.5,
            airfoil=vstab_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[
                -vstab_chord * vstab_taper_ratio ** 0.5,
                0,
                vstab_span
            ],
            chord=vstab_chord * vstab_taper_ratio ** 0.5,
            airfoil=vstab_airfoil,
        )
    ]
).translate([
    x_tail,
    0,
    0
])

### HStab
hstab_volume_coefficient = 0.47  # Potentially as low as 0.22
hstab_aspect_ratio = 4.5

hstab_area = hstab_volume_coefficient * wing.area() * wing.mean_aerodynamic_chord() / (x_tail - vstab_chord)
hstab_span = (hstab_area * hstab_aspect_ratio) ** 0.5
hstab_chord = (hstab_area / hstab_aspect_ratio) ** 0.5
hstab_airfoil = airfoils['ht14']

hstab = asb.Wing(
    name="HStab",
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[
                0,
                y,
                0,
            ],
            chord=hstab_chord,
            twist=4.00,
            airfoil=hstab_airfoil,
        )
        for y in [0, hstab_span / 2]
    ]
).translate([
    vstab.xsecs[-1].xyz_le[0],
    0,
    vstab_span,
])

### Fuselage
fuselage_length = (design_mass_TOGW / 4) ** (1 / 2)  # TODO optimize on this

# This controls the x-placement of the fuselage
x_fuse_nose = -0.35 * fuselage_length  # TODO link this to Cma stability calculation

fuse_wetted_area = 0.242552 * 2 * fuselage_length ** 2  # m^2, from CAD
fuse_xsec_area = 0.014823 * 2 * fuselage_length ** 2  # m^2, from CAD
fuse_volume = (8975383.23 * 1e-9) * 2 * fuselage_length ** 3  # m^3, from CAD

### Wing sponsons
wing_sponson_length = wing_root_chord * 0.75
wing_sponson_diameter = wing_sponson_length * 0.3
wing_sponson_wetted_area = wing_sponson_length * (np.pi * wing_sponson_diameter)
wing_sponson_volume = wing_sponson_length * (np.pi / 4) * wing_sponson_diameter ** 2

### Boom
boom_fwd_diameter = 0.03
boom_aft_diameter = 0.02
boom_length = x_tail
boom_wetted_area = boom_length * (
        np.pi *
        (boom_fwd_diameter + boom_aft_diameter) / 2
)

### Propulsion

n_propellers = 2

ideal_propeller_pressure_jump = 0.20 * cruise_op_point.dynamic_pressure()

propeller_diameter = np.minimum(
    (4 / np.pi * design_thrust_cruise_total / ideal_propeller_pressure_jump / n_propellers) ** 0.5,
    # Set by ideal propeller pressure jump
    fuselage_length * (10 * u.inch)  # Set by water clearance needs
)

propulsive_area_total = n_propellers * (np.pi / 4) * propeller_diameter ** 2

motor_y_placement = 0.0889 + 1 * propeller_diameter

### Assemble Airplane
airplane = asb.Airplane(
    wings=[
        wing,
        hstab,
        vstab
    ]
)

##### Section: Propulsion

propeller_tip_mach = 0.36  # From Dongjoon, 4/30/20
propeller_rads_per_sec = propeller_tip_mach * cruise_op_point.atmosphere.speed_of_sound() / (propeller_diameter / 2)
propeller_rpm = propeller_rads_per_sec * 30 / np.pi

motor_efficiency = 0.75  # RC-grade estimate

propeller_coefficient_of_performance = 0.90  # calibrated to QProp output with Dongjoon

cruise_power_propulsion = lib_prop_prop.propeller_shaft_power_from_thrust(
    thrust_force=design_thrust_cruise_total,
    area_propulsive=propulsive_area_total,
    airspeed=cruise_op_point.velocity,
    rho=cruise_op_point.atmosphere.density(),
    propeller_coefficient_of_performance=propeller_coefficient_of_performance
) / motor_efficiency
climb_power_propulsion = lib_prop_prop.propeller_shaft_power_from_thrust(
    thrust_force=design_thrust_climb_total,
    area_propulsive=propulsive_area_total,
    airspeed=cruise_op_point.velocity,
    rho=cruise_op_point.atmosphere.density(),
    propeller_coefficient_of_performance=propeller_coefficient_of_performance
) / motor_efficiency

# Motor thermal modeling
heat_motor = climb_power_propulsion * (1 - motor_efficiency) / n_propellers

propeller_max_torque = (climb_power_propulsion / n_propellers) / 0.7 / propeller_rads_per_sec

motor_kv = propeller_rpm / battery_voltage

##### Section: Solar Power Arrangement

### Panel arrangement

n_panels = 2 * (
        n_panels_spanwise_inboard * n_panels_chordwise_inboard +
        n_panels_spanwise_outboard * n_panels_chordwise_outboard
)

solar_area = n_panels * solar_area_per_panel

MPPT_efficiency = 1 / 1.04

power_in_at_panels = (
        0.5 * solar_area * lib_solar.solar_flux(
    latitude=latitude,
    day_of_year=day_of_year,
    time=time_after_solar_noon,
    altitude=cruise_op_point.atmosphere.altitude,
    panel_azimuth_angle=135 - 180,
    panel_tilt_angle=wing_dihedral_angle_deg,
) + 0.5 * solar_area * lib_solar.solar_flux(
    latitude=latitude,
    day_of_year=day_of_year,
    time=time_after_solar_noon,
    altitude=cruise_op_point.atmosphere.altitude,
    panel_azimuth_angle=135,
    panel_tilt_angle=wing_dihedral_angle_deg,
)
)

power_in_total = power_in_at_panels * MPPT_efficiency * solar_cell_efficiency

avionics_power = 8  # Watts

power_out_total_cruise = cruise_power_propulsion + avionics_power
power_out_total_climb = climb_power_propulsion + avionics_power

##### Section: Internal Geometry and Weights

structural_mass_markup = 1.2  # over Enigma F5J

mass_props = {}

### Lifting bodies
mass_wing = (  # Engima F5J
                    (0.440 + 0.460) *
                    (wing.area() / (0.264 * 3.624 * np.pi / 4)) ** 0.758 *
                    (design_mass_TOGW / 1.475) ** 0.49 *
                    (wing.aspect_ratio() / 18) ** 0.6
            ) * structural_mass_markup

mass_props['wing_center'] = asb.mass_properties_from_radius_of_gyration(
    mass=mass_wing * wing_y_break_fraction,
    x_cg=0,  # quarter-chord,
    radius_of_gyration_x=(wing_y_break_fraction * wing_span) / 12,
    radius_of_gyration_z=(wing_y_break_fraction * wing_span) / 12
)
mass_props['wing_tips'] = asb.mass_properties_from_radius_of_gyration(
    mass=mass_wing * (1 - wing_y_break_fraction),
    x_cg=0,  # quarter-chord,
    radius_of_gyration_x=(1 + wing_y_break_fraction) / 2 * (wing_span / 2),
    radius_of_gyration_z=(1 + wing_y_break_fraction) / 2 * (wing_span / 2),
)
mass_props['h_stab'] = asb.mass_properties_from_radius_of_gyration(
    mass=(
                 0.055 *
                 (
                         0.4 * (design_mass_TOGW / 1.475) ** 0.40 * (hstab_span / (0.670 / np.cosd(40))) ** 1.58 +
                         0.6 * (wing_root_chord / 0.264) * (hstab_span / (0.670 / np.cosd(40)))
                 )
         ) * structural_mass_markup,
    x_cg=hstab.xsecs[0].xyz_le[0] + hstab_chord / 2,
    z_cg=vstab.aerodynamic_center()[2],
    radius_of_gyration_x=hstab_span / 12,
    radius_of_gyration_y=hstab.xsecs[0].xyz_le[0] + hstab_chord / 2,
    radius_of_gyration_z=hstab.xsecs[0].xyz_le[0] + hstab_chord / 2
)
mass_props['v_stab'] = asb.mass_properties_from_radius_of_gyration(
    mass=(
                 0.055 *
                 (
                         0.3 * (design_mass_TOGW / 1.475) ** 0.40 * (vstab_span / (0.670 / np.cosd(40))) ** 1.58 +
                         0.7 * (wing_root_chord / 0.264) * (vstab_span / (0.670 / np.cosd(40)))
                 )
         ) * structural_mass_markup,
    x_cg=vstab.xsecs[0].xyz_le[0] + vstab_chord / 2,
    z_cg=vstab.aerodynamic_center()[2],
    radius_of_gyration_x=vstab_span / 12,
    radius_of_gyration_y=vstab.xsecs[0].xyz_le[0] + vstab_chord / 2,
    radius_of_gyration_z=vstab.xsecs[0].xyz_le[0] + vstab_chord / 2
)

### Other Structure
mass_props['boom'] = asb.mass_properties_from_radius_of_gyration(
    mass=(
                 0.235 *
                 (x_tail / 1.675) *
                 (design_mass_TOGW / 1.475) ** 0.49
         ) * structural_mass_markup,
    x_cg=x_tail / 2,
    radius_of_gyration_y=x_tail / 3,
    radius_of_gyration_z=x_tail / 3,
)

### Propulsion
mass_props['motors'] = asb.mass_properties_from_radius_of_gyration(
    mass=lib_prop_elec.mass_motor_electric(
        max_power=climb_power_propulsion / n_propellers,
        kv_rpm_volt=motor_kv,
        voltage=battery_voltage,
    ) * n_propellers,
    x_cg=-0.25 * wing_root_chord - 2 * u.inch,
    radius_of_gyration_x=motor_y_placement,
    radius_of_gyration_y=0.25 * wing_root_chord + 2 * u.inch,
    radius_of_gyration_z=motor_y_placement,
)

mass_props['motor_mounts'] = copy.copy(
    mass_props['motors']
) * 1  # similar to a quote from Raymer, modified to make sensible units, prop weight roughly subtracted

mass_props['propellers'] = asb.mass_properties_from_radius_of_gyration(
    mass=n_propellers * lib_prop_prop.mass_hpa_propeller(
        diameter=propeller_diameter,
        max_power=climb_power_propulsion / n_propellers,
        include_variable_pitch_mechanism=False
    ),
    x_cg=mass_props['motors'].x_cg - 1 * u.inch,
    radius_of_gyration_x=motor_y_placement,
    radius_of_gyration_z=motor_y_placement,
)

mass_props['ESCs'] = asb.mass_properties_from_radius_of_gyration(
    mass=lib_prop_elec.mass_ESC(max_power=climb_power_propulsion / n_propellers) * n_propellers,
    x_cg=0,
)

### Fuselage internals
mass_props['fuselage_skin'] = asb.mass_properties_from_radius_of_gyration(
    mass=(
                 fuse_wetted_area * (
                 np.minimum(2, 2 * fuselage_length) *  # plies
                 (4 * u.oz / u.yard ** 2)
         )
         ) * 1.5  # waterproofing markup
    ,
    x_cg=0.4 + x_fuse_nose,  # 0.4 to weight more towards thicker nose
    z_cg=-0.15 * fuselage_length
)
mass_props['fuselage_bulkheads'] = asb.mass_properties_from_radius_of_gyration(
    mass=(
            fuse_xsec_area * (  # These won't be anywhere near solid, but this is a rough estimate of material needed
            4 *  # n_bulkheads
            (1 / 8 * u.inch) * fuselage_length * (400)  # basswood
    )
    ),
    x_cg=(0.35 + x_fuse_nose) / 2,  # halfway between mid-fuse location and wing quarter-chord.
    z_cg=-0.15 * fuselage_length
)

mass_props['avionics'] = asb.mass_properties_from_radius_of_gyration(
    mass=0.060,  # RX, pixhawk mini
    x_cg=5 * u.inch + x_fuse_nose,
    z_cg=-0.15 * fuselage_length
)

mass_props['servos'] = asb.mass_properties_from_radius_of_gyration(
    mass=0.050,
    x_cg=x_tail * 0.5  # 2x tail, 2x wing
)

battery_capacity = (
        battery_reserve_time * (3 / 4) * cruise_power_propulsion +
        battery_reserve_time * (1 / 4) * climb_power_propulsion
)

battery_ampacity_amp_hours = battery_capacity / u.hour / battery_voltage

mass_props['battery'] = asb.mass_properties_from_radius_of_gyration(
    mass=battery_capacity / u.hour / battery_specific_energy_Wh_kg,
    x_cg=x_fuse_nose + 3 * u.inch,
    z_cg=-0.15 * fuselage_length
)

mass_props['MPPTs'] = asb.mass_properties_from_radius_of_gyration(
    mass=(n_panels / 36) * 0.080,  # Genasun GV-5
    x_cg=x_fuse_nose + 5 * u.inch,
    z_cg=-0.15 * fuselage_length
)

mass_props['solar_cells'] = asb.mass_properties_from_radius_of_gyration(
    mass=solar_area * rho_solar_cells,
    x_cg=0.25 * wing_root_chord,  # at the half-chord (remembering that datum is at quarter-chord)
)

mass_props['wiring'] = asb.mass_properties_from_radius_of_gyration(
    mass=(n_panels / 72) * 0.100,
    x_cg=-wing_root_chord / 4,  # at the leading edge
    z_cg=-0.05 * fuselage_length
)

### Other
mass_props['wing_sponsons'] = asb.mass_properties_from_radius_of_gyration(
    mass=2 * (  # 2 sponsons
            wing_sponson_wetted_area * (
            2 *  # plies
            (4 * u.oz / u.yard ** 2)
    )
    ) * 1.5  # waterproofing markup
    ,
    x_cg=0  # quarter-chord
)
mass_props['wing_sponson_mounts'] = copy.copy(
    mass_props['wing_sponsons']
) * 2  # a total guess

### Summation
mass_props_TOGW = asb.MassProperties(mass=0)
for k, v in mass_props.items():
    mass_props_TOGW = mass_props_TOGW + v

### Add glue weight
mass_props['glue_weight'] = mass_props_TOGW * 0.08
mass_props_TOGW += mass_props['glue_weight']

##### Section: Aerodynamics

ab_aero = asb.AeroBuildup(
    airplane=airplane,
    op_point=cruise_op_point,
    xyz_ref=mass_props_TOGW.xyz_cg
).run_with_stability_derivatives(
    alpha=True,
    beta=True,
    p=False,
    q=False,
    r=False
)

"""
Below are the results of a CFD study:
* Geometry: Fuselage scaled to 1 meter length
* Flow conditions:
    * 101325 Pa, 15 C air properties
    * 14 m/s freestream airspeed
    * 0 degrees angle of attack
"""
fuse_cfd_op_point = asb.OperatingPoint(velocity=14)
fuse_cfd_results = dict(
    CDA=0.171317 * 2 / fuse_cfd_op_point.dynamic_pressure() * fuselage_length ** 2,  # m^2
    CLA=0.0455037 * 2 / fuse_cfd_op_point.dynamic_pressure() * fuselage_length ** 2,  # m^2
    momentAc=0.0499564 * 2 / fuse_cfd_op_point.dynamic_pressure() * fuselage_length ** 3,  # m^3
)

fuse_flat_plate_results = dict(
    drag=(
            lib_aero.Cf_flat_plate(fuse_cfd_op_point.reynolds(reference_length=1)) *
            fuse_wetted_area *
            fuse_cfd_op_point.dynamic_pressure()
    )
)

aero = {}

aero['L'] = ab_aero['L'] + 0.
aero['D'] = ab_aero['D'] + (
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
aero['m_b'] = ab_aero['m_b'] + fuse_cfd_results['momentAc'] * cruise_op_point.dynamic_pressure()

aero['CL'] = aero['L'] / cruise_op_point.dynamic_pressure() / airplane.s_ref
aero['CD'] = aero['D'] / cruise_op_point.dynamic_pressure() / airplane.s_ref
aero['Cm'] = aero['m_b'] / cruise_op_point.dynamic_pressure() / airplane.s_ref / airplane.c_ref

aero['LD_ideal'] = aero['L'] / aero['D']

aero['Cma_fuse'] = (
        2 * (fuse_volume + 2 * wing_sponson_volume) / airplane.s_ref / airplane.c_ref
)
aero['Cma'] = ab_aero['Cma'] + aero['Cma_fuse']

aero['Cnb_fuse'] = (
        -2 * (fuse_volume + 2 * wing_sponson_volume) / airplane.s_ref / airplane.b_ref
)
aero['Cnb'] = ab_aero['Cnb'] + aero['Cnb_fuse']

##### Section: Compute other quantities
excess_power_cruise = power_in_total - power_out_total_cruise
breakeven_climb_rate = (
        excess_power_cruise * cruise_op_point.velocity /
        (cruise_power_propulsion * LD_cruise)
)

breakeven_climb_angle_deg = np.arcsind(breakeven_climb_rate / cruise_op_point.velocity)

##### Section: Constrain Problem
opti.subject_to([
    design_mass_TOGW > mass_props_TOGW.mass,  # mass closure
    LD_cruise < 0.75 * aero['LD_ideal'],  # aerodynamics closure, with knockdown to account for imperfect flying,
    aero['CL'] > 0,  # This keeps the velocity from getting confused
    aero['CL'] < 0.8,  # can't be riding the stall line too hard
    aero['L'] / 9.81 > design_mass_TOGW  # airplane must lift itself
])
opti.minimize(-breakeven_climb_rate)

if __name__ == '__main__':
    try:
        sol = opti.solve()
    except RuntimeError:
        sol = opti.debug
    s = lambda x: sol.value(x)

    airplane.substitute_solution(sol)
    cruise_op_point.substitute_solution(sol)
    # dyn.substitute_solution(sol)
    mass_props_TOGW.substitute_solution(sol)

    for v in mass_props.values():
        v.substitute_solution(sol)

    aero = {
        k: s(v)
        for k, v in aero.items() if not isinstance(v, list)
    }

    avl_aero = asb.AVL(
        airplane=airplane,
        op_point=cruise_op_point,
        xyz_ref=mass_props_TOGW.xyz_cg
    ).run()

    import matplotlib.pyplot as plt
    import aerosandbox.tools.pretty_plots as p
    from aerosandbox.tools.string_formatting import eng_string

    ##### Section: Printout
    print_title = lambda s: print(s.upper().join(["*" * 20] * 2))


    def fmt(x):
        return f"{s(x):.6g}"


    print_title("Outputs")
    for k, v in {
        "mass_TOGW"             : f"{fmt(mass_props_TOGW.mass)} kg ({fmt(mass_props_TOGW.mass / u.lbm)} lbm)",
        "L/D (actual)"          : fmt(LD_cruise),
        "Cruise Airspeed"       : f"{fmt(cruise_op_point.velocity)} m/s",
        "Cruise AoA"            : f"{fmt(cruise_op_point.alpha)} deg",
        "Cruise CL"             : fmt(aero['CL']),
        "breakeven_climb_rate"  : f"{fmt(breakeven_climb_rate)} m/s ({fmt(breakeven_climb_rate / (u.foot / u.minute))} ft/min)",
        "breakeven_climb_angle" : f"{fmt(breakeven_climb_angle_deg)} deg",
        "Cma"                   : fmt(aero['Cma']),
        "Cnb"                   : fmt(aero['Cnb']),
        "Cm"                    : fmt(aero['Cm']),
        "Wing Reynolds Number"  : eng_string(cruise_op_point.reynolds(wing.mean_aerodynamic_chord())),
        "AVL: Cma + Cma_fuse"   : avl_aero['Cma'] + aero['Cma_fuse'],
        "AVL: Cnb + Cnb_fuse"   : avl_aero['Cnb'] + aero['Cnb_fuse'],
        "AVL: Cm"               : avl_aero['Cm'],
        "AVL: Clb Cnr / Clr Cnb": avl_aero['Clb Cnr / Clr Cnb'] * avl_aero['Cnb'] / (
                avl_aero['Cnb'] + aero['Cnb_fuse']),
        "CG location"           : "(" + ", ".join([fmt(xyz) for xyz in mass_props_TOGW.xyz_cg]) + ") m",
        "Wing Span"             : f"{fmt(wing_span)} m ({fmt(wing_span / u.foot)} ft)",
        "Propeller Diameter"    : f"{fmt(propeller_diameter)} m ({fmt(propeller_diameter / u.inch)} in)",
        "S_ailerons / S_wing"   : fmt(aileron_area_fraction),
    }.items():
        print(f"{k.rjust(25)} = {v}")

    fmtpow = lambda x: fmt(x) + " W"

    print_title("Powers")
    for k, v in {
        "power_in_total"        : fmtpow(power_in_total),
        "power_out_total_cruise": fmtpow(power_out_total_cruise),
        "power_out_total_climb" : fmtpow(power_out_total_climb),
        "n_panels"              : int(np.round(s(n_panels))),
    }.items():
        print(f"{k.rjust(25)} = {v}")

    print_title("Mass props")
    for k, v in mass_props.items():
        print(f"{k.rjust(25)} = {v.mass:.3f} kg ({v.mass / u.oz:.2f} oz)")

    if make_plots:
        ##### Section: Geometry
        airplane.draw_three_view(show=False)
        p.show_plot(tight_layout=False, savefig="figures/three_view.png")

        ##### Section: Mass Budget
        fig, ax = plt.subplots(figsize=(12, 5), subplot_kw=dict(aspect="equal"), dpi=300)

        name_remaps = {
            **{
                k: k.replace("_", " ").title()
                for k in mass_props.keys()
            },
            "apu"                         : "APU",
            "wing"                        : "Wing Structure",
            "hstab"                       : "H-Stab Structure",
            "vstab"                       : "V-Stab Structure",
            "fuselage"                    : "Fuselage Structure",
            "payload_proportional_weights": "FAs, Food, Galleys, Lavatories,\nLuggage Hold, Doors, Lighting,\nAir Cond., Entertainment"
        }

        p.pie(
            values=[
                v.mass
                for v in mass_props.values()
            ],
            names=[
                n if n not in name_remaps.keys() else name_remaps[n]
                for n in mass_props.keys()
            ],
            center_text=f"$\\bf{{Mass\\ Budget}}$\nTOGW: {s(mass_props_TOGW.mass):.3f} kg",
            label_format=lambda name, value, percentage: f"{name}, {value:.3f} kg, {percentage:.1f}%",
            startangle=110,
            arm_length=30,
            arm_radius=20,
            y_max_labels=1.1
        )
        p.show_plot(savefig="figures/mass_budget.png")
