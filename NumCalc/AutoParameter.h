//
// File: AutoParameter.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov 11 22:15:16 2003
//

#ifndef _AUTOPARAMETER_H_
#define _AUTOPARAMETER_H_

#include "Parameter.h"

/**
 * @brief The AutoParameter class.
 *
 * This class overides the setValue() method of class Parameter so that no
 * Exception is thrown. This allows to perform optimization under constraint.
 */ 

class AutoParameter : public Parameter
{
	protected:
		ostream * _messageHandler;
		static double TINY;
	
	public: // Class constructors and destructors:
		/**
		 * @brief Build a new AutotParameter.
		 *
		 * @param name The parameter name.
		 * @param value The parameter value.
		 * @param constraint An optional pointer toward a constraint object.
		 */
		AutoParameter(const string & name = "", double value = 0, const Constraint * constraint = NULL);

		//Copy constructor:
		/**
		 * @brief Copy constructor.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter(const Parameter & param);
	
		/**
		 * @brief Copy constructor.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter(const AutoParameter & param);

		/**
		 * @brief Assignment operator.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter & operator=(const Parameter & param);
	
		//Destructor:
		virtual ~AutoParameter();
	
	public:
		
		Clonable * clone() const;
	
		/**
		 * @brief Set the value of this parameter.
		 *
		 * This method is overdefined so that no constraintException is thrown!
		 * When a Constraint is match, then we automatically apply a correct value instead.
		 * This correct value is the nearer limit reached by the value, or a value next to
		 * the limit if the limit is not reachable.
		 *
		 * This allow to perform optimisation under constraint whith algorithm that are not
		 * initially build for this.
		 *
		 * @param value the new parameter value.
		 */
		virtual void setValue(double value) throw (ConstraintException);
	
	public: //Specific method:
		
		/**
		 * @brief Set the message handler for this AutoParameter.
		 *
		 * The message handler keeps all messages that the paramter may send.
		 * The default handler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * @param mh The message handler to use.
		 */
		virtual void setMessageHandler(ostream * mh);
		
};


#endif	//_AUTOPARAMETER_H_
