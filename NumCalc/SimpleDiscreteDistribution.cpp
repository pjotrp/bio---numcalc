#include "SimpleDiscreteDistribution.h"

// From Utils:
#include <Utils/MapTools.h>

SimpleDiscreteDistribution::SimpleDiscreteDistribution(
	const map<double, double> & distribution
) {
	for(map<double, double>::const_iterator i = distribution.begin(); i != distribution.end(); i++)
		_distribution[i -> first] = i -> second;
}

SimpleDiscreteDistribution::~SimpleDiscreteDistribution() {}

Domain SimpleDiscreteDistribution::getDomain() const
{
	// Compute a new arbitray bounderi:
	vector<double> values = MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(_distribution);
	unsigned int n = values.size(); 
	vector<double> bounderi(n + 1);
	
	// Fill from 1 to n-1 with midpoints:
	for(unsigned int i = 1; i <= n - 1; i++)
		bounderi[i] = (values[i] - values[i - 1]) / 2.;
	
	// Fill 0 with the values[0] - (midpoint[0] - values[0]):
	bounderi[0] = 2 * values[0] - bounderi[1];
	
	// Fill n with values[n - 1] + (values[n - 1] - midpoint[n - 1]):
	bounderi[n] = 2 * values[n - 1] - bounderi[n - 1];
	
	// Build a domain and return it
	return Domain(bounderi, values);
}
