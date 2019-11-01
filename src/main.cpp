#include "simulation.h"
#include "random.h"

/*!
  \mainpage %Neuron network dynamics

  This is an implementation of the model of E.M. Izhikevich
  [(Simple Model of Spiking Neuron, IEE Trans. Neural Net., 2003)](https://www.izhikevich.org/publications/spikes.pdf).

  A typical command is 
  \verbatim
 ./NeuronNet -T IB:0.2,CH:0.15 -N 800 -d 100 -o test1000 -t 1000
  \endverbatim
  wich will create a \ref Network "network" of 800 \ref Neuron "neurons" consisting of 20% of *IB* neurons, 15% of *CH* neurons and *RS* neurons (see Neuron::NeuronTypes). 

  \ref Network::random_connect "Random connections" will be created so that the average neuron has 100 incoming connections. A simulation of a 1000 time-steps will be \ref Simulation::run "run" and the results will be saved in 3 \ref Simulation::output "output" files:
 - test1000: raster plots (table of 1000 columns and 800 lines with 1 when a neuron *n* fires at time *t*, 0 elsewhere),
 - test1000_traj: one representative time trajectory for each neuron type,
 - test1000_pars: all neuron parameters.

 */

RandomNumbers *_RNG;

int main(int argc, char **argv) {
    _RNG = new RandomNumbers;
    try {
        Simulation Sim(argc, argv);
        Sim.run();
    } catch(TCLAP::ArgException &e) {
        throw(TCLAP_ERROR("Error: " + e.error() + " " + e.argId()));
    } catch (SimulError &e) {
        std::cerr << e.what() << std::endl;
        return e.value();
    }
    if (_RNG) delete _RNG;
    return 0;
}
