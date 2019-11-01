#include "globals.h"
#include "neuron.h"

/*! \class Network
  A neuron network is a \ref neurons "set" of neurons and a \ref links "set" of directional links between them.

  Neurons are objects of class Neuron (they are identified by their index in the vector \ref neurons).
  A link is an ordered pair of indices in \ref neurons with an intensity value, 
  collected in a \ref linkmap. This is a [std::map](https://en.cppreference.com/w/cpp/container/map) 
  therefore only one connection can exist between two neurons. 

  To create a network, you need to \ref resize it and optionally \ref set_default_params for each neuron. 
  Then you can either call \ref add_link for each connection or generate a random network with \ref random_connect.

  The dynamics of the network proceeds by calling \ref step. 
  The state of the network can be printed to output streams with \ref print_params (to print all parameters of all neurons), \ref print_traj to print the full state of one neuron of each type. 
  The helper function \ref print_head will print a header line with the variable names for the the \ref print_traj lines.

  This class provides accessors to the following Neuron properties
  - \ref degree : returns the degree of a neuron (number of incoming connections) and valence (sum of weights of these links)
  - \ref neighbors : returns the indices of neurons with incoming links to a given neuron,
  - \ref potentials : returns the values of membrane potentials for all neurons,
  - \ref recoveries : returns the values of recovery variables for all neurons,
  - \ref neuron : returns a const reference to a given neuron object.

 */

typedef std::map<std::pair<size_t, size_t>, double> linkmap;

class Network {

public:
/*! 
  Resizes a network (grow or shrink). 
  \param n (size_t): new size of the network. If growing it will be filled with default excitatory and inhibitory neurons.
  \param _i (double): proportion of inhibitory neurons among the newly created.
 */
    void resize(const size_t&, double _i=_PROP_INHIB_);
/*! 
  Sets the neurons parameters.
  \param _types : a map between \ref Neuron::NeuronTypes and the number of neurons to set.
  \param _s : neurons will be set consecutively in \ref neurons starting from that index.
 */
    void set_default_params(const std::map<std::string, size_t>&, const size_t _s=0);
/*! 
  Sets the neurons parameters to non-default values.
  \param _types : a vector of \ref Neuron::NeuronTypes.
  \param _par : a vector of parameter values (\ref NeuronParams).
  \param _start : index of starting neuron in \ref neurons.
 */
    void set_types_params(const std::vector<std::string>&,
                          const std::vector<NeuronParams>&,
                          const size_t _s=0);
/*! 
  Sets the values of all membrane potentials.
  \param _poten : a vector of potential Neuron::_poten values.
  \param _s : index of starting neuron in \ref neurons.
 */
    void set_values(const std::vector<double>&, const size_t _s=0);
/*! 
  Creates a new link in \ref links.
  \param a (size_t): receiving neuron,
  \param b (size_t): sending neuron,
  \param str (double): link intensity (will be multiplied by -2 for inhibitory source).
  \return true if the link could be created.
 */
    bool add_link(const size_t&, const size_t&, double);
/*! 
  Creates random links in the network. Each neuron will expect to receive *n* connections (RandomNumbers::poisson with mean \p mean_deg) of intensity *s* (RandomNumbers::uniform_double with mean \p mean_streng).
  Sending neurons will be picked at random and since there can be only one connection from each neuron, the expected degree is not always achieved.
  \param mean_deg (double): mean value of Poisson distribution.
  \param mean_streng (double): mean value of the uniform distribution (with bounds 0 and 2*mean_streng).
  \return the number of links created.
 */
    size_t random_connect(const double&, const double &s=_STRENG_);
    size_t size() const {return neurons.size();}
/*! 
  Calculates the number and total intensity of connections to neuron \p n.
  \param n : the index of the receiving neuron.
  \return a pair {number of connections, sum of link intensities}.
 */
    std::pair<size_t, double> degree(const size_t&) const; //
    const Neuron& neuron(const size_t n) const {return neurons.at(n);}
/*! 
  Finds the list of neurons with incoming connections to \p n.
  \param n : the index of the receiving neuron.
  \return a vector of pairs {neuron index, link intensity}.
 */
    std::vector<std::pair<size_t, double> > neighbors(const size_t&) const;//
    std::vector<double> potentials() const;
    std::vector<double> recoveries() const;
/*! 
  Performs one time-step of the simulation.
  \param input : a vector of random values as thalamic input, one value for each neuron. The variance of these values corresponds to excitatory neurons.
  \return the indices of firing neurons.
 */
    std::set<size_t> step(const std::vector<double>&);//
    void print_params(std::ostream *_out=&std::cout);
    void print_traj(const int, const std::map<std::string, size_t>&, 
                    std::ostream *_out=&std::cout);
    void print_head(const std::map<std::string, size_t>&, 
                    std::ostream *_out=&std::cout);

private:
    std::vector<Neuron> neurons;
    linkmap links;

};
