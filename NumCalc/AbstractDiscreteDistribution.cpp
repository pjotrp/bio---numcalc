#include "AbstractDiscreteDistribution.h"

#include "VectorTools.h"
#include "RandomTools.h"

using namespace VectorFunctions;

/******************************************************************************/
	
unsigned int AbstractDiscreteDistribution::getNumberOfCategories() const {
	return _distribution.size();
}

/******************************************************************************/

double AbstractDiscreteDistribution::getCategory(unsigned int categoryIndex) const
{ 
	map<double, double>::const_iterator it = _distribution.begin();
	for(unsigned int i = 0; i < categoryIndex; i++) it++;
	return it -> first;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(unsigned int categoryIndex) const
{
	map<double, double>::const_iterator it = _distribution.begin();
	for(unsigned int i = 0; i < categoryIndex; i++) it++;
	return it -> second;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(double category) const 
{
	return _distribution.find(category) -> second;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getCategories() const
{
	Vdouble result(_distribution.size());
	unsigned int i = 0;
	for(map<double, double>::const_iterator it = _distribution.begin();
		it != _distribution.end();
		it++)
	{
		result[i] = it -> first;
		i++;
	}
	return result;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getProbabilities() const
{
	Vdouble result(_distribution.size());
	int i = 0;
	for(map<double, double>::const_iterator it = _distribution.begin();
		it != _distribution.end();
		it++) 
	{
		result[i] = it -> second;
		i++;
	}
	return result;
}

/******************************************************************************/

void AbstractDiscreteDistribution::set(double category, double probability) {
	_distribution[category] = probability;
}

/******************************************************************************/

void AbstractDiscreteDistribution::add(double category, double probability) {
	if(_distribution.find(category) == _distribution.end()) {
		//new category
		_distribution[category] = probability;
	} else {
		//existing category
		_distribution[category] += probability;
	}
}

/******************************************************************************/

double AbstractDiscreteDistribution::rand() const 
{
	double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	double cumprob = 0;
	for(map<double,double>::const_iterator i = _distribution.begin(); 
		i != _distribution.end();
		i++)
	{
		cumprob += i -> second;
		if(r <= cumprob) return i -> first;
	}
	// This line can't be reached:
	return -1.;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getInfCumulativeProbability(double category) const 
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	for(map<double, double>::const_iterator i = _distribution.begin();
		i != it;
		i++) prob += i -> second;
	return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getIInfCumulativeProbability(double category) const
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	if(it == _distribution.end()) return 0;
	for(map<double, double>::const_iterator i = ++it;
		i != _distribution.end();
		i++) prob += i -> second;
	return 1. - prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSupCumulativeProbability(double category) const 
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	if(it == _distribution.end()) return 0;
	for(map<double, double>::const_iterator i = ++it;
		i != _distribution.end();
		i++) prob += i -> second;
	return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSSupCumulativeProbability(double category) const
{
	double prob = 0;
	map<double, double>::const_iterator it = _distribution.find(category);
	for(map<double, double>::const_iterator i = _distribution.begin(); 
		i != it;
		i++) prob += i -> second;
	return 1. - prob;
}

/******************************************************************************/

void AbstractDiscreteDistribution::print(ostream & out) const
{
	for(map<double, double>::const_iterator i = _distribution.begin(); i != _distribution.end(); i++) {
		out << "Pr(" << (i -> first) << ") = " << (i -> second) << endl;
	}
}

/******************************************************************************/

ParameterList AbstractDiscreteDistribution::getParameters() const {
	return _parameters;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getParameter(const string & name) const
throw (ParameterNotFoundException)
{
	Parameter * p = _parameters.getParameter(name);
	if(p == NULL) throw ParameterNotFoundException("AbstractDiscreteDistribution::getParameter", name);
	return p -> getValue();
}

/******************************************************************************/

void AbstractDiscreteDistribution::setAllParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setAllParametersValues(params);
	fireParameterChanged();	
}

/******************************************************************************/

void AbstractDiscreteDistribution::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParameterValue(name, value);
	fireParameterChanged();
}

/******************************************************************************/

void AbstractDiscreteDistribution::setParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParametersValues(params);
	fireParameterChanged();	
}

/******************************************************************************/

void AbstractDiscreteDistribution::matchParametersValues(const ParameterList & params)
throw (ConstraintException)
{
	_parameters.matchParametersValues(params);
	fireParameterChanged();	
}

/******************************************************************************/
