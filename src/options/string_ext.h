#ifndef STRING_EXT_H
#define STRING_EXT_H

#include <string>
#include <cstring>
#include <sstream>
#include <vector>

#include "common.h"

#ifdef __BORLANDC__
#	include <ctype.h>
#	include <stdio.h>

#define str_transform1(first1, last1, d_first, fun) { \
	std::string::iterator _first1 = first1; \
	std::string::iterator _last1 = last1; \
	std::string::iterator _d_first = str.begin(); \
	while (_first1 != _last1) *_d_first++ = fun(*_first1++); \
}

#else
#	include <cstdio>
#include <algorithm>
#	define str_transform1 std::transform
#endif

#define SS_PRECISION 10

#ifdef STRING_EXT_USE_LOCALE
#define str2double str2double_l
#define num2str num2str_l
#else
#define str2double str2double_c
#define num2str num2str_c
#endif

// versions using localized decimal separator:
//double str2double_locale(std::string str);
double str2double_l(const std::string & str);
// versions using dot:
double str2double_c(const std::string & str);

template<class T>
std::string num2str_c(const T & value) {
	std::stringstream out;
	out.precision(SS_PRECISION);
	out << value;
	return out.str();
}

std::string num2str_l(const double & value);
std::string num2str_l(const int & value);
std::string num2str_l(const long int & value);


/*
template<class T>
std::string num2str_l(const T & value) {
	std::stringstream out;
	out.precision(SS_PRECISION);
	// locale support BEGIN
	const char* locstr = setlocale(LC_NUMERIC, NULL);
	std::locale loc = std::locale(locstr);
	delete locstr;
	out.imbue(loc);
	// locale support END

	out << value;
	return out.str();
}
*/

int str2int(const std::string & str);

void itrim(std::string & str, char const* sepSet = " \n\r\f");
void istrip_quotes(std::string & s);


//TODO: merge with itrim, istrip_quotes
void strtrim(std::string & str, const char * sepSet = " \n\t\r");
void strunquote(std::string & s);
// TODO: bool strip_quotes = , bool trimWs =
void strsplit(const std::string & in, std::string & left, std::string & right, const char c);





std::string trim(const std::string & str, char const* sepSet = " \n\r\f");

std::vector<std::string> str_tokenize(const std::string &source,
		const char *delimiter = " ", bool keepEmpty = false,
		bool trimWhitespace = true);

std::string string_printf(const char * fmt, ...)
              __attribute__ ((format (printf, 1, 2)));

extern inline std::string strtolower(std::string str) {
	str_transform1(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

extern inline std::string strtoupper(std::string str) {
	str_transform1(str.begin(), str.end(), str.begin(), ::toupper);
	return str;
}

extern inline char* str2charptr(const std::string & str) {
	size_t n = str.length();
	char* name = new char[n + 1];
	strncpy(name, str.c_str(), n);
	name[n] = 0;
	return name;
}

bool is_quoted_string(const std::string & s);
std::string strip_quotes(const std::string & s);


inline std::string strfill(int size, char ch) {
	std::string res(size, ch);
	return res;
}

std::string append_suffix_to_filename(const std::string & filename,
		const std::string & sfx);

void append_suffix_to_filename(char* out, const std::string & filename,
		const std::string & sfx);


#endif
