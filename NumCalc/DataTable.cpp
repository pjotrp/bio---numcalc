//
// File: DataTable.cpp
// Created by: Julien Dutheil
// Created on: Aug 2005
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "DataTable.h"
#include "VectorTools.h"

//From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/StringTokenizer.h>

using namespace bpp;

/******************************************************************************/

DataTable::DataTable(unsigned int nRow, unsigned int nCol) :
  _nRow(nRow), _nCol(nCol)
{
  _data.resize(nCol);
  for(unsigned int i = 0; i < nCol; i++)
    _data[i].resize(nRow);
  _colNames = NULL;
  _rowNames = NULL;
}
  
DataTable::DataTable(unsigned int nCol) :
  _nRow(0), _nCol(nCol)
{
  _data.resize(nCol);
  _colNames = NULL;
  _rowNames = NULL;
}
  
DataTable::DataTable(const vector<string> & colNames) throw (DuplicatedTableColumnNameException) :
  _nRow(0), _nCol(colNames.size())
{
  _data.resize(_nCol);
  _colNames = NULL;
  _rowNames = NULL;
  setColumnNames(colNames); //May throw an exception.  
}
  
DataTable::DataTable(const DataTable & table) :
  _nRow(table._nRow), _nCol(table._nCol), _data(table._data), _rowNames(NULL), _colNames(NULL)
{
  if(table._rowNames) _rowNames = new vector<string>(*table._rowNames);
  if(table._colNames) _colNames = new vector<string>(*table._colNames);
}

DataTable & DataTable::operator=(const DataTable & table)
{
  _nRow = table._nRow;
  _nCol = table._nCol;
  _data = table._data;
  if(_rowNames) delete _rowNames;
  if(_colNames) delete _colNames;
  _rowNames = NULL;
  _colNames = NULL;
  if(table._rowNames) _rowNames = new vector<string>(*table._rowNames);
  if(table._colNames) _colNames = new vector<string>(*table._colNames);
  return *this;
}

/******************************************************************************/

DataTable::~DataTable()
{
  if(_rowNames != NULL) delete _rowNames;
  if(_colNames != NULL) delete _colNames;
}

/******************************************************************************/
/*                             Cell access                                    */
/******************************************************************************/

string & DataTable::operator()(unsigned int rowIndex, unsigned int colIndex) throw (IndexOutOfBoundsException)
{
  if(colIndex >= _nCol) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, unsigned int).", colIndex, 0, _nCol - 1);
  if(rowIndex >= _data[colIndex].size()) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, unsigned int).", rowIndex, 0, _data[colIndex].size() - 1);
  return _data[colIndex][rowIndex];
}

const string & DataTable::operator()(unsigned int rowIndex, unsigned int colIndex) const throw (IndexOutOfBoundsException)
{
  if(colIndex >= _nCol) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, unsigned int).", colIndex, 0, _nCol - 1);
  if(rowIndex >= _data[colIndex].size()) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, unsigned int).", rowIndex, 0, _data[colIndex].size() - 1);
  return _data[colIndex][rowIndex];
}

/******************************************************************************/

string & DataTable::operator()(const string & rowName, const string & colName)
throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::operator(const string &, const string &).");
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::operator(const string &, const string &).");
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

const string & DataTable::operator()(const string & rowName, const string & colName) const
throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::operator(const string &, const string &).");
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::operator(const string &, const string &).");
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex) 
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}
    
/******************************************************************************/

string & DataTable::operator()(const string & rowName, unsigned int colIndex)
throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::operator(const string &, unsigned int).");
  if(colIndex >= _nCol) throw IndexOutOfBoundsException("DataTable::operator(const string &, unsigned int).", colIndex, 0, _nCol - 1);
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, unsigned int).", *ex.getElement());
  }
}

const string & DataTable::operator()(const string & rowName, unsigned int colIndex) const
throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::operator(const string &, unsigned int).");
  if(colIndex >= _nCol) throw IndexOutOfBoundsException("DataTable::operator(const string &, unsigned int).", colIndex, 0, _nCol - 1);
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, unsigned int).", *ex.getElement());
  }
}

/******************************************************************************/

string & DataTable::operator()(unsigned int rowIndex, const string & colName)
throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::operator(unsigned int, const string &).");
  try
  {
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    if(rowIndex >= _data[colIndex].size()) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, const string &).", rowIndex, 0, _data[colIndex].size() - 1);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

const string & DataTable::operator()(unsigned int rowIndex, const string & colName) const
throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::operator(unsigned int, const string &).");
  try
  {
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    if(rowIndex >= _data[colIndex].size()) throw IndexOutOfBoundsException("DataTable::operator(unsigned int, const string &).", rowIndex, 0, _data[colIndex].size() - 1);
    return (*this)(rowIndex, colIndex);
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

/******************************************************************************/
/*                             Work with names                                */
/******************************************************************************/

void DataTable::setRowNames(const vector<string> & rowNames)
throw (DimensionException, DuplicatedTableRowNameException)
{
  if(!VectorTools::isUnique(rowNames))
  {
    throw DuplicatedTableRowNameException("DataTable::setRowNames(...). Row names must be unique.");
  }
  if(rowNames.size() != _nRow) throw DimensionException("DataTable::setRowNames.", rowNames.size(), _nRow);
  else
  {
    if(_rowNames != NULL) delete _rowNames;
    _rowNames = new vector<string>(rowNames.begin(), rowNames.end());
  }
}

vector<string> DataTable::getRowNames() const throw (NoTableRowNamesException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::getRowNames().");
  return * _rowNames;
}

string DataTable::getRowName(unsigned int index) const throw (NoTableRowNamesException, IndexOutOfBoundsException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::getRowName(unsigned int).");
  if(index >= _nRow)    throw IndexOutOfBoundsException("DataTable::getRowName(unsigned int).", index, 0, _nRow - 1);
  return (* _rowNames)[index];
}

/******************************************************************************/

void DataTable::setColumnNames(const vector<string> & colNames)
throw (DimensionException, DuplicatedTableColumnNameException)
{
  if(!VectorTools::isUnique(colNames)) throw DuplicatedTableColumnNameException("DataTable::setColumnNames(...). Column names must be unique.");
  if(colNames.size() != _nCol) throw DimensionException("DataTable::setColumnNames.", colNames.size(), _nCol);
  else
  {
    if(_colNames != NULL) delete _colNames;
    _colNames = new vector<string>(colNames.begin(), colNames.end());
  }
}

vector<string> DataTable::getColumnNames() const throw (NoTableColumnNamesException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::getColumnNames().");
  return *_colNames;
}

string DataTable::getColumnName(unsigned int index) const throw (NoTableColumnNamesException, IndexOutOfBoundsException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::getColumnName(unsigned int).");
  if(index >= _nCol)    throw IndexOutOfBoundsException("DataTable::getColumnName(unsigned int).", index, 0, _nCol - 1);
  return (* _colNames)[index];
}

/******************************************************************************/
/*                               Work on columns                              */
/******************************************************************************/

vector<string> & DataTable::getColumn(unsigned int index)
  throw (IndexOutOfBoundsException)
{
  if(index >= _nCol) throw IndexOutOfBoundsException("DataTable::getColumn(unsigned int).", index, 0, _nCol - 1);
  return _data[index];
}

const vector<string> & DataTable::getColumn(unsigned int index) const
  throw (IndexOutOfBoundsException)
{
  if(index >= _nCol) throw IndexOutOfBoundsException("DataTable::getColumn(unsigned int).", index, 0, _nCol - 1);
  return _data[index];
}  

vector<string> & DataTable::getColumn(const string & colName)
  throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::getColumn(const string &).");
  try
  {
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    return _data[colIndex];
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableColumnNameNotFoundException("DataTable::getColumn(const string &).", colName);
  }
}

const vector<string> & DataTable::getColumn(const string & colName) const
  throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::getColumn(const string &).");
  try
  {
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    return _data[colIndex];
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableColumnNameNotFoundException("DataTable::getColumn(const string &).", colName);
  }
}

bool DataTable::hasColumn(const string & colName) const
{ 
  if(_colNames == NULL) return false;
  for(unsigned int i = 0; i < _colNames->size(); i++)
  {
    if((* _colNames)[i] == colName) return true;
  }
  return false;
}

void DataTable::deleteColumn(unsigned int index)
  throw (IndexOutOfBoundsException)
{
  if(index >= _nCol) throw IndexOutOfBoundsException("DataTable::deleteColumn(unsigned int).", index, 0, _nCol - 1);
  _data.erase(_data.begin() + index);
  if(_colNames != NULL) _colNames->erase(_colNames->begin()+index);
  _nCol--;
}

void DataTable::deleteColumn(const string & colName)
  throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
{
  if(_colNames == NULL) throw NoTableColumnNamesException("DataTable::deleteColumn(const string &).");
  try
  {
    unsigned int colIndex = VectorTools::which(*_colNames, colName);
    _data.erase(_data.begin() + colIndex);
    _colNames->erase(_colNames->begin() + colIndex);
    _nCol--;
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableColumnNameNotFoundException("DataTable::deleteColumn(const string &).", colName);
  }
}

void DataTable::addColumn(const vector<string> & newColumn)
  throw (DimensionException, TableColumnNamesException)
{
  if(_colNames != NULL) throw TableColumnNamesException("DataTable::addColumn. Table has column names.");
  if(newColumn.size() != _nRow) throw DimensionException("DataTable::addColumn.", newColumn.size(), _nRow);
  for(unsigned int i = 0; i < _nRow; i++)
  {
    _data[i].push_back(newColumn[i]);
  }
  _nCol++;
}

void DataTable::addColumn(const string & colName, const vector<string> & newColumn)
  throw (DimensionException, NoTableColumnNamesException, DuplicatedTableColumnNameException)
{
  if(_colNames == NULL)
  {
    if(_nCol == 0) _colNames = new vector<string>();
    else throw NoTableColumnNamesException("DataTable::addColumn. Table has column names.");
  }
  if(newColumn.size() != _nRow) throw DimensionException("DataTable::addColumn.", newColumn.size(), _nRow);
  if(_nCol > 0 && find(_colNames->begin(), _colNames->end(), colName) != _colNames->end())
    throw DuplicatedTableColumnNameException("DataTable::addColumn(const string &, const vector<string> &). Column names must be unique.");
  _colNames->push_back(colName);
  _data.push_back(newColumn);
  _nCol++;
}

/******************************************************************************/
/*                               Work on rows                                 */
/******************************************************************************/

vector<string> DataTable::getRow(unsigned int index) const
  throw (IndexOutOfBoundsException)
{
  if(index >= _nRow) throw IndexOutOfBoundsException("DataTable::getRow(unsigned int).", index, 0, _nRow - 1);
  vector<string> row;
  for (unsigned int i = 0 ; i < _nCol ; i++) {
    row.push_back(_data[i][index]);
  }
  return row;
}

vector<string> DataTable::getRow(const string & rowName) const
  throw (NoTableRowNamesException, TableRowNameNotFoundException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::getRow(const string &).");
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    vector<string> row;
    for (unsigned int i = 0 ; i < _nCol ; i++)
      row.push_back(_data[i][rowIndex]);
    return row;
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableRowNameNotFoundException("DataTable::getRow(const string &).", rowName);
  }
}

bool DataTable::hasRow(const string & rowName) const
{ 
  if(_rowNames == NULL) return false;
  for(unsigned int i = 0; i < _rowNames->size(); i++) {
    if((* _rowNames)[i] == rowName) return true;
  }
  return false;
}

void DataTable::deleteRow(unsigned int index)
  throw (IndexOutOfBoundsException)
{
  for(unsigned int j = 0; j < _nCol; j++) {
    vector<string> * column = & _data[j];
    if(index >= column->size()) throw IndexOutOfBoundsException("DataTable::deleteRow(unsigned int).", index, 0, column->size() - 1);
    column->erase(column->begin()+index);
  }
  if(_rowNames != NULL) _rowNames->erase(_rowNames->begin()+index);
  _nRow--;
}

void DataTable::deleteRow(const string & rowName)
  throw (NoTableRowNamesException, TableRowNameNotFoundException)
{
  if(_rowNames == NULL) throw NoTableRowNamesException("DataTable::deleteRow(const string &).");
  try
  {
    unsigned int rowIndex = VectorTools::which(*_rowNames, rowName);
    for(unsigned int j = 0; j < _nCol; j++)
    {
      vector<string> * column = & _data[j];
      column->erase(column->begin()+rowIndex);
    }
    _rowNames->erase(_rowNames->begin()+rowIndex);
    _nRow--;
  }
  catch(ElementNotFoundException<string> & ex)
  {
    throw TableRowNameNotFoundException("DataTable::deleteRow(const string &).", rowName);
  }
}

void DataTable::addRow(const vector<string> & newRow)
  throw (DimensionException, TableRowNamesException)
{
  if(_rowNames != NULL) throw TableRowNamesException("DataTable::addRow. Table has row names.");
  if(newRow.size() != _nCol) throw DimensionException("DataTable::addRow.", newRow.size(), _nCol);
  for(unsigned int j = 0; j < _nCol; j++) {
    _data[j].push_back(newRow[j]);
  }
  _nRow++;
}

void DataTable::addRow(const string & rowName, const vector<string> & newRow)
  throw (DimensionException, NoTableRowNamesException, DuplicatedTableRowNameException)
{
  if(_rowNames == NULL) {
    if(_nRow == 0) _rowNames = new vector<string>();
    else throw NoTableRowNamesException("DataTable::addRow. Table has row names.");
  }
  if(newRow.size() != _nCol) throw DimensionException("DataTable::addRow.", newRow.size(), _nCol);
  if(_nRow > 0 && find(_rowNames->begin(), _rowNames->end(), rowName) != _rowNames->end())
      throw DuplicatedTableRowNameException("DataTable::addRow(const string &, const vector<string> &). Row names must be unique.");
  _rowNames->push_back(rowName);
  for(unsigned int j = 0; j < _nCol; j++) {
    _data[j].push_back(newRow[j]);
  }
  _nRow++;
}

/******************************************************************************/
/*                               Read from a CSV file                         */
/******************************************************************************/

DataTable * DataTable::read(istream & in, const string & sep, bool header, int rowNames)
  throw (DimensionException, IndexOutOfBoundsException, DuplicatedTableRowNameException)
{
  string firstLine  = FileTools::getNextLine(in);
  StringTokenizer st1(firstLine, sep);
  vector<string> row1(st1.getTokens().begin(), st1.getTokens().end());
  string secondLine = FileTools::getNextLine(in);
  StringTokenizer st2(secondLine, sep);
  vector<string> row2(st2.getTokens().begin(), st2.getTokens().end());
  unsigned int nCol = row1.size();
  bool hasRowNames;
  DataTable * dt;
  if(row1.size() == row2.size())
  {
    dt = new DataTable(nCol);
    if(header)
    { //Use first line as header.
      dt->setColumnNames(row1);
    }
    else
    {
      dt->addRow(row1);
    }
    dt->addRow(row2);
    hasRowNames = false;
  }
  else if(row1.size() == row2.size() - 1)
  {
    dt = new DataTable(nCol);
    dt->setColumnNames(row1);
    string rowName = *row2.begin();
    dt->addRow(rowName, vector<string>(row2.begin()+1, row2.end()));
    hasRowNames = true;
  }
  else throw DimensionException("DataTable::read(...). Row 2 has not the correct number of columns.", row2.size(), nCol);

  // Now read each line:
  string line = FileTools::getNextLine(in);
  while(!TextTools::isEmpty(line))
  {
    StringTokenizer st(line, sep);
    if(hasRowNames)
    {
      string rowName = *st.getTokens().begin();
      vector<string> row(st.getTokens().begin()+1, st.getTokens().end());
      dt->addRow(rowName, row);
    }
    else
    {
      vector<string> row(st.getTokens().begin(), st.getTokens().end());
      dt->addRow(row);
    }
    line = FileTools::getNextLine(in);
  }

  // Row names:
  if(rowNames > -1)
  {
    if((unsigned int)rowNames >= nCol) throw IndexOutOfBoundsException("DataTable::read(...). Invalid column specified for row names.", rowNames, 0, nCol-1);
    vector<string> col = dt->getColumn((unsigned int)rowNames);
    dt->setRowNames(col);
    dt->deleteColumn(rowNames);
  }
  
  return(dt);
}

/******************************************************************************/

void DataTable::write(const DataTable & data, ostream & out, const string & sep)
{
  unsigned int n = data.getNumberOfColumns();
  if(n == 0) return;
  if(data.hasColumnNames())
  { //Write header
    vector<string> colNames = data.getColumnNames();
    out << colNames[0];
    for(unsigned int i = 1; i < n; i++) {
      out << sep << colNames[i];
    }
    out << endl;
  }
  //Now write each row:
  for(unsigned int i = 0; i < data.getNumberOfRows(); i++)
  {
    if(data.hasRowNames()) out << data.getRowName(i) << sep;
    out << data(i, 0);
    for(unsigned int j = 1; j < n; j++)
    {
      out << sep << data(i, j);
    }
    out << endl;
  }
}

/******************************************************************************/

