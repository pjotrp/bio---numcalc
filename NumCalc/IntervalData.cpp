//
// File: IntervalData.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 19:59:29 2004
//

#include "IntervalData.h"

#include <cmath>

IntervalData::IntervalData(const Domain & domain, const string & name):
	_domain(domain), _name(name)
{
	reset();
}

Domain IntervalData::getDomain() const { return _domain; }

double IntervalData::getDomainValue(double x) const { return _domain.getNearestValue(x); }

void IntervalData::setName(const string & name) { _name = name; }

string IntervalData::getName() const { return _name; }

unsigned int IntervalData::getFreq(unsigned int i) const { return _freqs[i]; }

double IntervalData::getDensity(unsigned int i) const { return (double)_freqs[i] / (double)_n; }

vector<unsigned int> IntervalData::getFrequencies() const { return _freqs; }

Vdouble IntervalData::getDensities() const
{
	Vdouble densities(_freqs.size());
	for(unsigned int i = 0; i < _freqs.size(); i++) densities[i] = getDensity(i);
	return densities;
}

void IntervalData::addValue(double value) throw (OutOfRangeException)
{
	//_data.push_back(value);
	_n++;
	_sum += value;
	_sumsquare += value * value;
	if(value < _min) _min = value;
	if(value > _max) _max = value;
	int index = _domain.getIndex(value);
	_freqs[index]++;
}

unsigned int IntervalData::getSize() const { return _n; }

double IntervalData::getMinValue() const { return _min; }

double IntervalData::getMaxValue() const { return _max; } 
	
double IntervalData::getMean() const { return _sum / _n; } 
		
double IntervalData::getSD() const { return (_n / (_n - 1)) * getSDP(); } 
	
double IntervalData::getSDP() const { return _sumsquare / _n - _sum * _sum / (_n * _n); }

void IntervalData::reset()
{
	_freqs.resize(_domain.getSize());
	_sum = 0;
	_sumsquare = 0;
	_min = -log(0.);
	_max = log(0.);
	_n = 0;
}

void IntervalData::print(ostream & out) const
{
	out << "midpoint\tlowerB\tupperB\tfreq\tdensity" << endl;
	for(unsigned int i = 0; i < _domain.getSize(); i++) {
		out << _domain.getValue(i);
		out << "\t";
		out << _domain.getBound(i);
		out << "\t";
		out << _domain.getBound(i + 1);
		out << "\t";
		out << _freqs[i];
		out << "\t";
		out << getDensity(i);
		out << endl;
	}
}
