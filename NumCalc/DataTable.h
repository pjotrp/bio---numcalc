/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _DataFrame_H_
#define _DataFrame_H_

#include "VectorTools.h"

// From Utils:
#include <Utils/Exceptions.h>
#include <Utils/TextTools.h>

// From the STL:
#include <string>
#include <map>
using namespace std;

class TableNameNotFoundException: public Exception
{
	protected:
		string _name;
		
	public:
		TableNameNotFoundException(const string & text, const string & name) :
			Exception("NameNotFoundException: "+text), _name(name) {}
		~TableNameNotFoundException() throw() {}

	public:
		string getName() const { return _name; }		
};

class NoTableRowNamesException: public Exception
{
	public:
		NoTableRowNamesException(const string & text) :
			Exception("NoTableRowNamesException: "+text) {}
		~NoTableRowNamesException() throw() {}
};

class NoTableColumnNamesException: public Exception
{
	public:
		NoTableColumnNamesException(const string & text) :
			Exception("NoTableColumnNamesException: "+text) {}
		~NoTableColumnNamesException() throw() {}
};

class TableRowNamesException: public Exception
{
	public:
		TableRowNamesException(const string & text) :
			Exception("TableRowNamesException: "+text) {}
		~TableRowNamesException() throw() {}
};

class TableColumnNamesException: public Exception
{
	public:
		TableColumnNamesException(const string & text) :
			Exception("TableColumnNamesException: "+text) {}
		~TableColumnNamesException() throw() {}
};

class DuplicatedTableRowNameException: public Exception
{
	public:
		DuplicatedTableRowNameException(const string & text) :
			Exception("DuplicatedTableRowNameException: "+text) {}
		~DuplicatedTableRowNameException() throw() {}
};

class DuplicatedTableColumnNameException: public Exception
{
	public:
		DuplicatedTableColumnNameException(const string & text) :
			Exception("DuplicatedTableColumnNameException: "+text) {}
		~DuplicatedTableColumnNameException() throw() {}
};




/**
 * @brief This class corresponds to a 'dataset', <i>i.e.</i> a table with data by rows and variable by columns.
 *
 * Data are stored as string, by column.
 * A DataTable object is hence similar to a ColMatrix<string>.object.
 */
class DataTable {
	
	protected:
		unsigned int _nRow, _nCol;
		vector< vector<string> > _data;
		vector<string> * _rowNames;
		vector<string> * _colNames;

	public:
		
		/**
		 * @brief Builds a new void DataTable object with nRow rows and nCol columns.
		 *
		 * @param nRow The number of rows of the DataTable.
		 * @param nCol The number of columns of the DataTable.
		 */
		DataTable(unsigned int nRow, unsigned int nCol);
		
		/**
		 * @brief Builds a new void DataTable object with nCol columns.
		 *
		 * @param nCol The number of columns of the DataTable.
		 */
		DataTable(unsigned int nCol);

		/**
		 * @brief Builds a new void DataTable object with named columns.
		 *
		 * @param colNames The names of the columns of the DataTable.
		 * @throw DuplicatedTableColumnNameException Ifcolnames contains identical names.
		 */
		DataTable(const vector<string> & colNames) throw (DuplicatedTableColumnNameException);

		~DataTable();

	public:
		      string & operator()(unsigned int rowIndex, unsigned int colIndex)       throw (IndexOutOfBoundsException);
		const string & operator()(unsigned int rowIndex, unsigned int colIndex) const throw (IndexOutOfBoundsException);
		      string & operator()(const string & rowName, const string & colName)
						throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException, DimensionException);
		const string & operator()(const string & rowName, const string & colName) const
						throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException, DimensionException);
		
		unsigned int getNumberOfRows() const { return _nRow; }
		unsigned int getNumberOfColumns() const { return _nCol; }
		void setRowNames(const vector<string> & rowNames) throw (DimensionException, DuplicatedTableRowNameException);
		vector<string> getRowNames() const throw (NoTableRowNamesException);
		string getRowName(unsigned int index) const throw (NoTableRowNamesException, IndexOutOfBoundsException);
		bool hasRowNames() const { return _rowNames!= NULL; }
		void setColumnNames(const vector<string> & colNames) throw (DimensionException, DuplicatedTableColumnNameException);
		vector<string> getColumnNames() const throw (NoTableColumnNamesException);
		string getColumnName(unsigned int index) const throw (NoTableColumnNamesException, IndexOutOfBoundsException);
		bool hasColumnNames() const { return _colNames!= NULL; }

    /**
		 * @name Work on columns.
		 *
		 * @{
		 */
		      vector<string> & getColumn(unsigned int index)       throw (IndexOutOfBoundsException);
		const vector<string> & getColumn(unsigned int index) const throw (IndexOutOfBoundsException);
		void deleteColumn(unsigned int index) throw (IndexOutOfBoundsException);
		void addColumn(const vector<string> & newColumn) throw (DimensionException, TableColumnNamesException);
		void addColumn(const string & colName, const vector<string> & newColumn) throw (DimensionException, NoTableColumnNamesException, DuplicatedTableColumnNameException);
		/** @} */
    
		/**
		 * @name Work on rows.
		 *
		 * @{
		 */
		void deleteRow(unsigned int index) throw (IndexOutOfBoundsException);
		void addRow(const vector<string> & newRow) throw (DimensionException, TableRowNamesException);
		void addRow(const string & rowName, const vector<string> & newRow) throw (DimensionException, NoTableRowNamesException, DuplicatedTableRowNameException);
		/** @} */

	public:
		/**
		 * @brief Read a table form a CSV file.
		 *
		 * The number of rows is given by the second line in the file.
		 * By default, if the first line as one column less than the second one,
		 * the first line is taken as column names, and the first column as row names.
		 * Otherwise, no column names and no row names are specified, unless
		 * explicitely precised by the user.
		 * 
		 * @param in       The input stream.
		 * @param sep      The column delimiter.
		 * @param header   Tell if the first line must be used as column names, otherwise use default.
		 * @param rowNames Use a column as rowNames. If positive, use the specified column to compute rownames, otherwise use default;
		 * @return         A pointer toward a new DataTable object.
		 */
		static DataTable * read(istream & in, const string & sep = "\t", bool header=true, int rowNames=-1)
			throw (DimensionException, IndexOutOfBoundsException);

		static void write(const DataTable & data, ostream & out, const string & sep = "\t");
};

#endif //_DataFrame_H_
