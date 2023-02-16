from design_opt import *
import aerosandbox.numpy as np

if __name__ == '__main__':

    npsi, npci = np.meshgrid(
        np.arange(1, 21),
        np.arange(1, 9)
    )

    npso = npsi * (3 / 2)
    npco = np.minimum(npci - 1, npci * 2 / 3)

    ws = (
            2 * (npsi + npso) * panel_spacing +
            wing_extra_span
    )

    sols = opti.solve_sweep(
        {
            n_panels_spanwise_inboard  : npsi,
            n_panels_chordwise_inboard : npci,
            n_panels_spanwise_outboard : npso,
            n_panels_chordwise_outboard: npco,
        },
        update_initial_guesses_between_solves=True,
        solve_kwargs=dict(
            max_iter=50,
        )
    )


    def val(x):
        return np.vectorize(
            lambda sol: np.nan if sol is None else sol.value(x),
            cache=True
        )(sols)


    import matplotlib.pyplot as plt
    import aerosandbox.tools.pretty_plots as p


    def finalize_plot(title=""):
        plt.annotate(
            text="Seaway-Mini",
            xy=(3.84807, 3),
            xytext=(4.5, 3.5),
            xycoords="data",
            arrowprops={
                "color"     : "k",
                "width"     : 0.25,
                "headwidth" : 4,
                "headlength": 6,
                "shrink"    : 0.1
            }
        )
        plt.plot([3.84807], [3], ".k")
        # plt.ylim(bottom=0)
        plt.xlim(right=8)
        p.show_plot(
            title,
            "Wing Span [m]",
            "# of Chordwise Panels in Center",
            tight_layout=False
        )


    fig, ax = plt.subplots()
    p.contour(
        ws,
        npci,
        val(breakeven_climb_angle_deg),
        colorbar_label="Breakeven Climb Gradient [deg]",
        linelabels_format=lambda x: f"${x:.1f}^\\circ$",
        levels=np.arange(0, 8)
    )
    finalize_plot("Net-Zero-Energy Climb Gradient\nOn April 1, Boston, 2 p.m.")

    # fig, ax = plt.subplots()
    # p.contour(
    #     ws,
    #     npci,
    #     val(breakeven_climb_rate),
    #     colorbar_label="Breakeven Climb Rate [m/s]",
    #     linelabels_format=lambda x: f"{x:.1f} m/s",
    #     levels=np.arange(0, 1.1, 0.1)
    # )
    # finalize_plot("Net-Zero-Energy Climb Rate\nOn April 1, Boston, 2 p.m.")

    fig, ax = plt.subplots()
    p.contour(
        ws,
        npci,
        val(excess_power_cruise / mass_props_TOGW.mass),
        colorbar_label="Excess Specific Energy [W/kg]",
        linelabels_format=lambda x: f"{x:.1f} W/kg",
        # levels=np.arange(0, 1.1, 0.1)
    )
    finalize_plot("Excess Specific Energy At Cruise\nOn April 1, Boston, 2 p.m.")

    fig, ax = plt.subplots()
    p.contour(
        ws,
        npci,
        val(design_mass_TOGW),
        colorbar_label="TOGW [kg]",
        linelabels_format=lambda x: f"{x:.1f} kg",
        z_log_scale=True,
    )
    finalize_plot("Expected Takeoff Gross Weight")

    mass_struct = sum([
        mass_props[m].mass for m in [
            "wing_center",
            "wing_tips",
            "h_stab",
            "v_stab",
            "boom",
            "fuselage_skin",
            "fuselage_bulkheads",
            "wing_sponsons",
            "wing_sponson_mounts",
            "glue_weight"
        ]
    ])

    fig, ax = plt.subplots()
    p.contour(
        ws,
        npci,
        val(mass_struct/design_mass_TOGW),
        colorbar_label="Structural Mass Fraction [%]",
        linelabels_format=lambda x: f"{x*100:.0f}%",
        z_log_scale=True,
    )
    finalize_plot("Structural Mass Fraction")

    fig, ax = plt.subplots()
    p.contour(
        ws,
        npci,
        val(LD_cruise),
        colorbar_label="$L/D$ (actual)",
        linelabels_format=lambda x: f"{x:.1f}",
        levels=np.arange(4, 22)
    )
    finalize_plot("Aerodynamic Efficiency $L/D$")
