/* JAM first attempt at a toy two pool model with HMM kinetics in
 Rust.  This generalised model should be expanadable to any number of
 pools and interactions.  I'll use rk4 integration algorithm only
 because it is what I have used historically and it worked!  These
 biological systems, comprised of Henri-Michaelis-Menten (HMM) kinetic
 equations are usually not stiff.  This model follows the structure of
the accompanying diagram called "Two Pool Model.pdf */


// First produced by Pablo Zamora, Sieglord, and JAM at Norwich UK on
// 2026_01_26

//Last updated on 2026_02_01

use russell_lab::NumVector;
use russell_ode::{Method, OdeSolver, Params, System};
use gnuplot::*;
use std::thread::sleep;
use std::time::Duration;

// --- Model Parameters ---

// Number of state variables (e.g., amount of metabolite in Pool A and Pool B).
const N_STATES: usize = 2;

// Capacities or initial parameters for the system.
// I[0], I[1]: Volumes/capacities for Pool A and B.
// I[4]: Constant external input to Pool A.

//SA=20.0, SB=25.0,QA0= 6.0,QB0=9.0, FOA=3.0
const I:[f64;5] = [20.0, 25.0,6.0,9.0,3.0];

/* Declaration of kinetic constants for the HMM equations.  A V12
  variable refers to a VMax for the equation describing a flux from
  pool "1" to pool "2".  A K12 variable refers to the "affinity
  constant" for the equation describing a flux from pool "1" to pool
  "2". */

// Flux constants used in the Michaelis-Menten-style equations.
// VAB = 18.0, VBA = 13.0, VBO = 8.0, KAB = 0.32,
// KBA = 0.36, KBO = 0.31
const C: [f64; 6] = [18.0,13.0,8.0,0.32,0.36,0.31];

// --- Non-State Variable Extraction ---

/// This struct holds "non-state" variables (results like concentrations and fluxes)
/// that are calculated during the ODE integration but are not part of the state vector.
#[derive(Debug, Clone, Copy)]
pub struct AuxiliaryResults {
    pub con_a: f64,
    pub con_b: f64,
    pub fab: f64,
    pub fba: f64,
    pub fbo: f64,
}

impl AuxiliaryResults {
    // Initializer for the auxiliary results.
    fn new() -> Self {
        AuxiliaryResults {
            con_a: I[3]/I[1],
            con_b: I[4]/I[2],
            fab: 0.0,
            fba: 0.0,
            fbo: 0.0,
        }
    }
}

fn main() {
    // 1. Define the system of differential equations
    // The mutable closure parameter 'results' allows the solver to update
    // auxiliary variables during each calculation stage.
    let system = System::new(N_STATES, |dydt, _t, y, results: &mut AuxiliaryResults| {
        // --- Calculate Concentrations ---
        results.con_a = y[0] / I[0];
        results.con_b = y[1] / I[1];

        // --- Calculate Fluxes (Mechanistic equations) ---
        results.fab = C[0] / (1.0 + (C[3] / results.con_a));
        results.fba = C[1] / (1.0 + (C[4] / results.con_b));
        results.fbo = C[2] / (1.0 + (C[5] / results.con_b));

        // --- Specify the ODEs ---
        dydt[0] = I[4] + results.fba - results.fab;
        dydt[1] = results.fab - results.fba - results.fbo;

        Ok(())
    });

    // 2. Configure the solver (Advance Runge-Kutta method)
    //   let params = Params::new(Method::DoPri8);
    // use Runge-Kutta 4th order let params =
    let params = Params::new(Method::Rk4);
    let mut solver = OdeSolver::new(params, system).expect("Solver
    initialization failed");

      // 3. Set Initial Conditions
    // start time
    let mut t = 0.0;
    // integration interval
    let dt = 0.1;
    // Initial metabolite amounts in pools
    let mut y = NumVector::from(&[I[4], I[3]]);

    let mut results = AuxiliaryResults::new();

    let mut fg = Figure::new();
    let num_steps = 100;
    let mut state_trace = Vec::with_capacity(num_steps + 1);
    let mut results_trace = Vec::with_capacity(num_steps + 1);
    let mut t_trace = Vec::with_capacity(num_steps + 1);
    t_trace.push(t);
    state_trace.push(y.clone());
    results_trace.push(results);

    // 4. Time-stepping Loop
    println!("Time, PoolA, PoolB, ConA, ConB, Fab, Fba, Fbo");
    //Number of steps 
    for _ in 0..num_steps {
        // Advance the simulation.
        // The 'results' struct is passed mutably so it captures the values
        // calculated inside the system function at the end of the step.
        solver
	//How would I know what to fill in the () for .solve?
            .solve(&mut y, t, t + dt, None, &mut results)
            .expect("Solver failed");
        t += dt;

        t_trace.push(t);
        state_trace.push(y.clone());
        results_trace.push(results);

        // Plot the trace.
        fg.clear_axes();
        fg.axes2d()
            .set_pos_grid(3, 1, 0)
            .set_x_range(Fix(0.), Fix(num_steps as f64 * dt))
            .lines_points(&t_trace, state_trace.iter().map(|y| y[0]), &[Caption("A")])
            .lines_points(&t_trace, state_trace.iter().map(|y| y[1]), &[Caption("B")])
            .lines_points(&t_trace, state_trace.iter().map(|y| y[0] + y[1]), &[Caption("Total")]);
        fg.axes2d()
            .set_pos_grid(3, 1, 1)
            .set_x_range(Fix(0.), Fix(num_steps as f64 * dt))
            .lines_points(&t_trace, results_trace.iter().map(|r| r.con_a), &[Caption("Con A")])
            .lines_points(&t_trace, results_trace.iter().map(|r| r.con_b), &[Caption("Con B")]);
        fg.axes2d()
            .set_pos_grid(3, 1, 2)
            .set_x_range(Fix(0.), Fix(num_steps as f64 * dt))
            .lines_points(&t_trace, results_trace.iter().map(|r| r.fab), &[Caption("Flux AB")])
            .lines_points(&t_trace, results_trace.iter().map(|r| r.fba), &[Caption("Flux BA")])
            .lines_points(&t_trace, results_trace.iter().map(|r| r.fbo), &[Caption("Flux BO")]);
        fg.show_and_keep_running().unwrap();
        sleep(Duration::from_millis(50));

        // Print the State variables (Pools) and Non-State variables (Cons/Fluxes)
        println!(
            "{:.2}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}",
	    // time,QA, QB, ConA,ConB,Fab,Fba,Fbo
            t, y[0], y[1], results.con_a, results.con_b, results.fab, results.fba, results.fbo
        );
    }
}
