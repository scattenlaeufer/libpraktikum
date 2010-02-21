#ifndef DATA
#define DATA

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

using namespace std;

/** \class Data
 * \brief Stores data of physical experiments.
 *
 */
class Data
{
	public:
		/** \brief Creates an empty data object
 		 * \param[in] length The length of the data
		 * \param[in] cols The number of columns
		 */
		Data(const unsigned int length, const unsigned int cols);

		/** \brief Creates a data object from a cassy file.
 		 * \param[in] filename The path to the cassy file.
		 */
		Data(const string &filename);

		/** Searches a column for the minimal value. The result is cached in the object, so this function has a minimal overhead after the first call.
 		 * \param[in] n the column
		 * \return the minimal value in the n column
		 */
		double getMin(const unsigned int n);

		/** Searches a column for the maximal value. The result is cached in the object, so this function has a minimal overhead after the first call.

 		 * \param[in] n the column
		 * \return the maximal value in the n column
		 */
		double getMax(const unsigned int n);

		/** \brief Get the n-th column of the data
 		 */
		inline double *getCol(const unsigned int n) {
			if (n < cols)
				return data[n];
			else
				return NULL;
		}

		/** \brief Get the values between two limits.
		 *
		 * \param[in] col The column which is searched
		 * \param[in] first The first limit
		 * \param[in] second The second limit
		 * \param[out] count The number of the found values
		 * \returns pointer to the array with the value between the two limits.
		 */
		double *getValuesBetween(const unsigned int col, const double first, const double second, int &count);

		/** \brief Get the name of the n-th column
		 */
		inline string getName(const unsigned int n) {
			if (names != NULL && n < cols)
				return names[n];
			else
				return NULL;
		}

		/** \brief Get the unit of the n-th column
 		 */
		inline string getUnit(const unsigned int n) {
			if (units != NULL && n < cols)
				return units[n];
			else
				return NULL;
		}

		/** \brief Get the symbol of the n-th column
 		 */
		inline string getSymbol(const unsigned int n) {
			if (symbols != NULL && n < cols)
				return symbols[n];
			else
				return NULL;
		}

		/** Get the length of the data. This is the length of each column.
 		 */
		inline unsigned int length() {
			return _length;
		}

	protected:
		/// The number of datasets
		unsigned int _length;
		/// The number of columns
		unsigned int cols;
		bool hasHeader;
		int posData;

		/// The table for the data
		double **data;
		/// The names of the columns
		string *names;
		/// The symbols, which are mostly used in formulas of the columns
		string *symbols;
		/// The units of the columns
		string *units;

		/// The cached minimal values for each column
		double *min;
		/// The cached maximal values for each column
		double *max;

		 /** \brief Scans Labfiles generated by cassy lab, to gather information about the amount of data
		 * 
		 * \note The number of cols and the number of data is stored in the labfile in the last line before the data block. \n
		 * This information is not taken into account to be able to process files with header removed or with lines removed.
		 * 
		 * \param[in] filename c++ string with filename of the labfile
		 * \param[out] &length returns numbers of found dataset
		 * \param[out] &cols returns numbers of data-cols found (determines the size of names, symbol unit and data)
		 * \param[out] &hasHeader returns weather a header was found or not (if not, names, symbol and unit will remain empty)
		 * \param[out] &posData position of the data in the data stream
		 * \returns 0 success \n
		 * 		1 Could not open input file \n
		 * 		2 No data found in input file
		 */
		int scanLab(string filename, unsigned int &length, unsigned int &cols, bool &hasHeader, int &posData);


		/**
		 * \brief Reads Header of labfile in already allocated arrays
		 * 
		 * \warning This might or might not work for a particular labfile
		 * 
		 * \param[in] filename c++ string with filename of the labfile
		 * \param[out] cols number of data-cols (can be determined with scanLab() )
		 * \param[out] names pointer to array with n_cols string elements, the names of the data columns will be put here
		 * \param[out] symbol pointer to array with n_cols string elements, an array with the symbols of the data coluns used in cassy lab will be put here
		 * \param[out] unit pointer to array with n_cols string elements, array with units for the data in each column will be given back here
		 * 
		 * \returns 0 Data read successfully \n
		 * 		1 Could not open input file
		 */

		int readLabHeader(string filename, unsigned int cols, string* names, string* symbol, string* unit);

		/**
		 * \brief Reads data from labfiles generated by cassy lab
		 * 
		 * \param[in] filename c++ string with filename of the labfile
		 * \param[out] &length numbers of found dataset (might be changed if a line turns out to not contain data)
		 * \param[out] cols numbers of data-cols 
		 * \param[out] data pointer to pointer to double, two dimensional array for data has to be allocated (data[n_cols][n_datasets])
		 * \param[in] posData position of the data in the data stream
		 * 
		 * \returns 0 Data read successfully \n
		 * 		1 Could not open input file
		 */

		int readLabData(string filename, unsigned int &length, unsigned int cols, double** data, int posData);

};



#endif
