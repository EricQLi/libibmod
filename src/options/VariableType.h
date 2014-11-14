/*
 * VariableType.h
 *
 *  Created on: 18 Feb 2013
 *      Author: s12kb2
 */

#ifndef VARIABLETYPE_H_
#define VARIABLETYPE_H_

//#include "defines.h"

#include <string>
#include <sstream>
//#include <math.h>

#ifdef __GNUC__
#include <cmath>
#else
#include <algorithm>
#define fmax std::max
#define fmin std::min
#endif

#include "string_ext.h"

#ifdef __INTEL_COMPILER
#define _return_type_const
#else
#define _return_type_const const
#endif


#define SEP_TYPE const std::string &

class NumericRange {
private:
	bool _integer;
	double lower, upper;
	static std::string separator;
public:
	NumericRange(double lower_, double upper_) {
		set(lower_, upper_);
		_integer = false;
	}

	NumericRange(int lower_, int upper_) {
		set((double) lower_, (double) upper_);
		_integer = true;
	}

	NumericRange(void) {
		set(0.0, 1.0);
		_integer = false;
	}

	NumericRange(std::string s, SEP_TYPE sep = separator) {
		set_from_string(s, sep);
		_integer = false;
	}
	NumericRange(const char* s, SEP_TYPE sep = separator) {
		set_from_string(s, sep);
		_integer = false;
	}

	void set(double lower_, double upper_) {
		lower = fmin(lower_, upper_);
		upper = fmax(lower_, upper_);
	}

	std::string to_string_c(SEP_TYPE sep = separator) const {
		std::stringstream ss;
		ss << lower << sep << upper;
		return ss.str();
	}

	std::string to_string(SEP_TYPE sep = separator) const {
		return num2str(lower) + sep + num2str(upper);
	}

	std::string to_string_l(SEP_TYPE sep = separator) const {
		return num2str_l(lower) + sep + num2str_l(upper);
	}

	void set_from_string(const char* s, SEP_TYPE sep = separator) {
		set_from_string(std::string(s), sep);
	}

	void set_from_string(const std::string & s, SEP_TYPE sep = separator) {
		int pos = s.find_first_of(sep);
		set(str2double(s.substr(0, pos)), str2double(s.substr(pos + 1)));
	}

	bool encloses(const double x) const {
		return x <= upper && x >= lower;
	}

	double fit_into(const double x) {
		if(x > upper) return upper;
		if(x < lower) return lower;
		return x;
	}

	static void set_separator(const std::string &sep) {
		separator = sep;
	}

	const std::string get_separator(void) const {
		return separator;
	}

	double get_lower(void) const {
		return lower;
	}

	double get_upper(void) const {
		return upper;
	}
	double diff(void) const {
		return upper - lower;
	}

	double mean(void) const {
		return (upper + lower) / 2.;
	}
};

#undef SEP_TYPE

class Factor {
private:
	std::vector<std::string> labels;
	size_t nlevels, value;
public:
	Factor(const std::string & strval, const std::string & labelstr, const char* delimiter);
	Factor(const int numval, const std::string & labelstr, const char* delimiter);

	Factor(void);

	~Factor();
	bool set(unsigned int v);
	bool set(const std::string & s);
	bool setLevels(const std::string & labelstr, const char* delimiter);

	std::string getLabel(void) const;
	int getLabel(char* out);

	std::string getLevel(size_t n) const;
	int getLevel(size_t n, char* out);

	int getNLevels(void) const {
		return nlevels;
	}

	int getInt(void) const;
};




class VariableType {
public:
	enum VAL_TYPE {
		VT_UNDEFINED = -1,
		VT_BOOLEAN = 1,
		VT_INTEGER = 2,
		VT_REAL = 3,
		VT_STRING = 4,
		VT_RANGE = 5,
		VT_FACTOR = 6
	};

private:
	VAL_TYPE type;
	long int intv;
	double doublev;
	std::string stringv;
	NumericRange rangev;
	Factor factorv;

public:
	VariableType(void) :
			type(VT_UNDEFINED), intv(0), doublev(0.0), stringv(""), rangev(), factorv() {
	}

	template<class T>
	VariableType(T val) :
			type(VT_UNDEFINED), intv(0), doublev(0.0), stringv(""), rangev(), factorv() {
		type = set(val);
	}

	std::string to_string(void) const;
	std::string to_string(bool use_locale) const;

	int get_int(void) const;
	double get_real(void) const;
	std::string get_str(void) const;
	bool get_bool(void) const;
	NumericRange get_range(void) const;
	Factor get_factor(void) const;

	VAL_TYPE get_type(void) const {
		return type;
	}

	void set_from_str(const std::string & val);
	void set_from_str(const char* val);

	VAL_TYPE set(const std::string & val);
	VAL_TYPE set(const char* val);
	VAL_TYPE set(const int val);
	VAL_TYPE set(const long int val);
	VAL_TYPE set(const double val);
	VAL_TYPE set(const bool val);
	VAL_TYPE set(const NumericRange & val);
	VAL_TYPE set(const Factor & val);

};

class VariableTypeOption: public VariableType {
private:
	VariableType default_data;
	std::string name;
	char short_name;
	unsigned int flags; // TODO: file/command line only
public:

	typedef enum {
		OPT_DEFAULT = -1,
		OPT_FILE = 1 /*1*/,
		OPT_CMDLINE = 1 << 1 /*2*/,
		OPT_FLAG = 1 << 2 /*4*/
	} OptionFlag;



	static int _counter; // XXX: move to private
	std::string description;
	const VariableType* get_default_value(void) const;
	bool validate(void) const;
	bool has_short_name(void) const {
		return (short_name != '\0');
	}

	const std::string get_name(void) const {
		return name;
	}
	_return_type_const char get_short_name(void) const {
		return short_name;
	}

	void _set_name(void) {
		unsigned int bar = name.find_first_of("|");
		if (bar != std::string::npos) {
			if (name.length() > bar)
				short_name = name[bar + 1];
			name = name.substr(0, bar);
		}
	}

	template<class T>
	VariableTypeOption(const T value, std::string _name,
			std::string _description) :
			VariableType(value), default_data(value), name(_name), short_name(
					'\0'), 
					flags(OPT_DEFAULT), 
					description(_description) {
		default_data.set(value);
		_set_name();
		_counter++;
	}

	VariableTypeOption(void) :
			VariableType(), default_data(), name("undefined"), short_name('\0'),
			 flags(OPT_DEFAULT), description("<undefined>") {
	}


	bool is_numeric(void) const {
		VAL_TYPE t = get_type();
		return t == VT_STRING || t == VT_BOOLEAN;
	}


};

#endif /* VARIABLETYPE_H_ */
