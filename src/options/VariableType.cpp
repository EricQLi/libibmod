/*
 * VariableType.cpp
 *
 *  Created on: 18 Feb 2013
 *      Author: s12kb2
 *      TODO: Handling of NAN, [+-]Inf
 *      TODO: Numeric levels for enum types
 */



#include "string_ext.h"
#include "VariableType.h"

#include <iostream>

#include "../defines.h"

std::string NumericRange::separator = "~";
int VariableTypeOption::_counter = 0;

Factor::Factor(const std::string & strval, const std::string & labelstr,
		const char* delimiter) : nlevels(0), value(0) {
	setLevels(labelstr, delimiter);
	if (!set(strval)) {
		printf("Error: in Factor, '%s' does not exist in set \n", strval.c_str());
		exit(1);
	}
}

Factor::Factor(const int numval, const std::string & labelstr,
		const char* delimiter) : nlevels(0), value(0) {
	setLevels(labelstr, delimiter);
	if (!set(numval)) {
		printf("Error: in Factor, numeric value (%d) is larger than number of levels (%u) \n", numval, (unsigned int) nlevels);
		exit(1);
	}
}



Factor::Factor(void) :
	nlevels(0), value(0) {
}


Factor::~Factor() {
}

//void Factor::_setLevels(const std::string & labelstr, const char* delimiter) {
//	labels = str_tokenize(labelstr, delimiter, false, true);
//	nlevels = labels.size();
//}

bool Factor::setLevels(const std::string & labelstr, const char* delimiter) {
	std::string strval;
	if(nlevels != 0) strval = labels[value];
	labels = str_tokenize(labelstr, delimiter, false, true);
	nlevels = labels.size();
	//printf("New factor: %d %s \n", (int) nlevels, labelstr.c_str());
	if(strval.empty()) return true;
	return set(strval);
}


bool Factor::set(unsigned int v) {
	if (v < nlevels) {
		value = v;
		return true;
	}
	return false;
}


bool Factor::set(const std::string & s) {
	value = std::find(labels.begin(), labels.end(), s) - labels.begin();
	return value < nlevels;
}

int Factor::getLevel(size_t n, char* out) {
	if (n < nlevels) {
		int len = labels[n].size();
		strcpy(out, labels[n].c_str());
		return len;
	}
	return 0;
}

std::string Factor::getLevel(size_t n) const {
	if (n < nlevels)
		return labels[n];
	return "";
}


std::string Factor::getLabel(void) const {
	return getLevel(value);
}

int Factor::getLabel(char* out) {
	return getLevel(value, out);
}


int Factor::getInt(void) const {
	if(nlevels == 0) return -1;
	return value;
}


VariableType::VAL_TYPE VariableType::set(const bool val) {
	intv = (long int) val;
	doublev = (double) val;
	stringv = val ? "True" : "False";
	return VT_BOOLEAN;
}
VariableType::VAL_TYPE VariableType::set(const int val) {
	intv = (long int) val;
	doublev = (double) val;
	stringv = num2str(intv);
	return VT_INTEGER;
}
VariableType::VAL_TYPE VariableType::set(const long int val) {
	intv = val;
	doublev = (double) val;
	stringv = num2str(intv);
	return VT_INTEGER;
}
VariableType::VAL_TYPE VariableType::set(const double val) {
	doublev = val;
	intv = (long int) val;
	stringv = num2str(doublev);
	return VT_REAL;
}

VariableType::VAL_TYPE VariableType::set(const char* val) {
	stringv = (std::string) val;
	doublev = str2double(stringv);
	intv = (long int) doublev;
	return VT_STRING;
}

VariableType::VAL_TYPE VariableType::set(const std::string & val) {
	stringv = (std::string) val;
	doublev = str2double(stringv);
	intv = (long int) doublev;
	return VT_STRING;
}
VariableType::VAL_TYPE VariableType::set(const NumericRange & val) {
	rangev = val;
	stringv = rangev.to_string();
	doublev = rangev.get_lower();
	intv = (long int) doublev;
	return VT_RANGE;
}

VariableType::VAL_TYPE VariableType::set(const Factor & val) {
	factorv = val;
	intv = factorv.getInt();
	stringv = factorv.getLabel();
	doublev = intv;
	return VT_FACTOR;
}




int VariableType::get_int(void) const {
	return intv;
}
double VariableType::get_real(void) const {
	return doublev;
}
std::string VariableType::get_str(void) const {
	return stringv;
}
bool VariableType::get_bool(void) const {
	return intv != 0;
}
NumericRange VariableType::get_range(void) const {
	return rangev;
}
Factor VariableType::get_factor(void) const {
	return factorv;
}


std::string VariableType::to_string(void) const {
	switch (type) {
	case VT_BOOLEAN:
	case VT_INTEGER:
	case VT_REAL:
	case VT_STRING:
	case VT_RANGE:
	case VT_FACTOR:
		return stringv;
		break;
	default:
		break;
	}
	return "<undefined>";
}

std::string VariableType::to_string(bool use_locale) const {
    switch(type) {
    case VT_REAL:
        return use_locale? num2str_l(doublev) : num2str_c(doublev);
    case VT_RANGE:
        return use_locale? rangev.to_string_l() : rangev.to_string_c();
    default:
        return to_string();
    }
}


void VariableType::set_from_str(const std::string & val) {
	std::string str = strip_quotes(trim(val));

	bool res;
	switch (type) {
	case VT_BOOLEAN:
		str = strtolower(trim(str, " \n\r\f\t"));
		res = ((str2double(str) != 0.0) || str == "true" || str == "yes"
				|| str == "y" || str == "1" || str == "t");
		set(res);
		break;
	case VT_INTEGER:
		set(str2int(str));
		break;
	case VT_REAL:
		set(str2double(str));
		break;
	case VT_STRING:
		set(val);
		break;
	case VT_RANGE:
		set(NumericRange(str));
		break;
	case VT_FACTOR:
//		XPRINT(str);
		factorv.set(str);
		stringv = factorv.getLabel();
		intv = factorv.getInt();
		doublev = (double) intv;
		break;
	default:
		break;
	}
}

void VariableType::set_from_str(const char* val) {
	set_from_str((const std::string) val);
}

const VariableType* VariableTypeOption::get_default_value(void) const {
//VariableType VariableTypeOption::get_default_value(void) {
	return &default_data;
}


