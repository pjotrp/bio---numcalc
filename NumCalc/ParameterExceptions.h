//
// File: ParameterExceptions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov  3 18:05:36 2003
//

#ifndef _PARAMETEREXCEPTIONS_H_
#define _PARAMETEREXCEPTIONS_H_

// From Utils:
#include <Utils/Exceptions.h>

class Parameter;

/******************************************************************************
 *                           Parameters exceptions:                           *
 ******************************************************************************/

/**
 * @brief The parameter exception base class.
 *
 * @see Exception
 */
class ParameterException : public Exception {

	protected:
		const Parameter * parameter;
			
	public:	// Class constructors and destructor:
		/**
		 * @brief Build a new ParameterException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param param A const pointer toward the parameter that threw the exception.
		 */	
		ParameterException(const char * text, const Parameter * param);
		/**
		 * @brief Build a new ParameterException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param param A const pointer toward the parameter that threw the exception.
		 */	
		ParameterException(const string & text, const Parameter * param);
	
		// Class destructor
		~ParameterException() throw ();
	public:
		/**
		 * @brief Get the parameter that threw the exception.
		 *
		 * @return A const pointer toward the parameter.
		 */
		virtual const Parameter * getParameter() const;
};

/**
 * @brief Exception thrown when a value is specified out of a constraint.
 */
class ConstraintException : public ParameterException {
	
	protected:
		double badValue;
	
	public: // Class constructor and destructor:
		/**
		 * @brief Build a new ConstraintException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    A const pointer toward the parameter that threw the exception.
		 * @param badValue The value that doesn't match the constraint.
		 */	
		ConstraintException(const char *   text, const Parameter * param, double badValue);
		/**
		 * @brief Build a new ConstraintException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    A const pointer toward the parameter that threw the exception.
		 * @param badValue The value that doesn't match the constraint.
		 */	
		ConstraintException(const string & text, const Parameter * param, double badValue);

		~ConstraintException() throw ();
	
	public:
		/**
		 * @brief Get the value that doesn't match the constraint.
		 *
		 * @return The faulty value.
		 */
		virtual double getBadValue() const;
};

/*******************************************************************************/

/**
 * @brief Exception thrown when a parameter is not found,
 * for instance in a ParameterList.
 */
class ParameterNotFoundException : public Exception {

	protected:
		const string parameter;
			
	public:	// Class constructors and destructor:
		/**
		 * @brief Build a new ParameterNotFoundException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    The name of the parameter not found.
		 */	
		ParameterNotFoundException(const char *   text, const string & param = "");
	
		/**
		 * @brief Build a new ParameterNotFoundException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    The name of the parameter not found.
		 */	
		ParameterNotFoundException(const string & text, const string & param = "");
	
		// Class destructor
		~ParameterNotFoundException() throw ();
	public:
		/**
		 * @brief Get the name of the parameter not found.
		 *
		 * @return The parameter name.
		 */
		virtual string getParameter() const;
};

#endif	//_PARAMETEREXCEPTIONS_H_
