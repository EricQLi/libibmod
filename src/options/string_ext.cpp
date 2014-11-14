#include <cstdarg>

#ifdef STRING_EXT_USE_LOCALE
#include <locale>
#endif

#include <cstdlib>

#include "string_ext.h"

std::string num2str_l(const double & value) {
	char res[64];
	sprintf(res, "%.10g", value);
	return std::string(res);
}

std::string num2str_l(const int & value) {
	char res[64];
	sprintf(res, "%d", value);
	return std::string(res);
}

std::string num2str_l(const long int & value) {
	char res[64];
	sprintf(res, "%ld", value);
	return std::string(res);
}

double str2double_c(const std::string & str) {
	if (((int) str.find_first_of("1234567890")) == -1)
		return 0;
	std::istringstream ss(str);
	//ss.getloc().name()
	//ss.imbue(std::locale("C"));
	double ret;
	ss >> ret;
	return ret;
}

double str2double_l(const std::string & str) {
	if (((int) str.find_first_of("1234567890")) == -1)
		return 0;
	return std::strtod(str.c_str(), NULL);
}


int str2int(const std::string & str) {
	if (((int) str.find_first_of("1234567890")) == -1)
		return 0;
#ifdef WITH_SSTREAM
	std::istringstream ss(str);
	int ret;
	ss >> ret;
	return ret;
#else
	return std::strtol(str.c_str(), NULL, 0);
#endif
}




/*double str2double_locale(std::string str) {
 if (((int) str.find_first_of("1234567890")) == -1)
 return 0.0;
 std::string cdec = ".";
 lconv* locale = localeconv();
 std::string sys_decimal_point = locale->decimal_point;
 if (sys_decimal_point != cdec) {
 unsigned int pos = str.find(sys_decimal_point);
 if (pos != std::string::npos)
 str.replace(pos, sys_decimal_point.length(), cdec);
 }
 std::stringstream ss;
 double ret;
 ss << str;
 ss >> ret;
 return ret;
 }*/

void itrim(std::string & str, const char * sepSet) {
	std::string::size_type const first = str.find_first_not_of(sepSet);
	if (first == std::string::npos)
		str = std::string();
	else
		str = str.substr(first, str.find_last_not_of(sepSet) - first + 1);
}


std::string trim(const std::string & str, const char * sepSet) {
	std::string::size_type const first = str.find_first_not_of(sepSet);
	return (first == std::string::npos) ?
			std::string() :
			str.substr(first, str.find_last_not_of(sepSet) - first + 1);
}


//Source: http://stackoverflow.com/a/10051869 (by Mahmoud Al-Qudsi, modified by KB)
std::vector<std::string> str_tokenize(const std::string & source,
		const char *delimiter, bool keepEmpty, bool trimWhitespace) {
	std::vector<std::string> results;

	std::string item;
	size_t prev = 0;
	size_t next = 0;

	while ((next = source.find_first_of(delimiter, prev)) != std::string::npos) {
		if (keepEmpty || (next - prev != 0)) {
			item = source.substr(prev, next - prev);
			if (trimWhitespace) itrim(item);
			if (keepEmpty || !item.empty())
				results.push_back(item);
		}
		prev = next + 1;
	}

	if (prev < source.size()) {
		item = source.substr(prev);
		if (trimWhitespace)	itrim(item);
		if (keepEmpty || !item.empty())
			results.push_back(item);
	}
	return results;
}



std::string string_printf(const char * fmt, ...) {
	std::string res;
	char* tmp_buf = NULL;
	long tmp_buf_size = 128;
	for (;;) {
		tmp_buf = (char *) malloc(tmp_buf_size);
		if (tmp_buf == NULL)
			break;
		va_list args;
		va_start(args, fmt);
		const long buf_size_needed = vsnprintf(tmp_buf, tmp_buf_size, fmt,
				args);
		va_end(args);
		if ((buf_size_needed) >= 0 && (buf_size_needed < tmp_buf_size)) {
			tmp_buf_size = buf_size_needed;
			break;
		}
		free(tmp_buf);
		if (buf_size_needed < 0) {
			tmp_buf_size *= 2;
		} else {
			tmp_buf_size = buf_size_needed + 1;
		}
	}
	if (tmp_buf != NULL) {
		res.append(tmp_buf, tmp_buf_size);
		free(tmp_buf);
	}
	return res;
}

/**
 * Tests whether the string represents a (single or double) quoted text.
 * @param s	the string
 * @return true if the string is quoted
 */
bool is_quoted_string(const std::string & s) {
	return (s.length() > 1) && ((s[0] == '"') || (s[0] == '\'')) && (s[0] == s[s.size() - 1]);
}

/**
 * Strips quotes (either single or double) from the given string
 * @param s the string
 * @return the string s with quotes removed
 */
std::string strip_quotes(const std::string & s) {
	return (is_quoted_string(s)) ? trim(s, "\"\'") : s;
}


void istrip_quotes(std::string & s) {
	if (is_quoted_string(s))
		itrim(s, "\"\'");
}


void strtrim(std::string & str, const char * sepSet) {
	std::string::size_type const first = str.find_first_not_of(sepSet);
	if (first == std::string::npos)
		str = std::string();
	else
		str = str.substr(first, str.find_last_not_of(sepSet) - first + 1);
}

void strunquote(std::string & s) {
	if ((s.length() > 1) && ((s[0] == '"') || (s[0] == '\'')) && (s[0] == s[s.size() - 1]))
		strtrim(s, "\"\'");
}


// TODO: bool strip_quotes = , bool trimWs =
void strsplit(const std::string & in, std::string& left, std::string& right, const char c) {
	size_t pos = in.find(c);
	if (pos == std::string::npos) {
		left = in;
		strtrim(left);
		right = "";
	} else if (pos <= 1) {
		left = "";
		right = in.substr(pos + 1, std::string::npos);
		strtrim(right);
		strunquote(right);
	} else {
		left = in.substr(0, pos);
		strtrim(left);
		right = in.substr(pos + 1, std::string::npos);
		strtrim(right);
		strunquote(right);
	}
}


/**
 * Appends a suffix to file base name
 * @param filename a string giving the file name or full path
 * @param sfx the suffix to append to base name
 * @return modified 'filename'
 */
std::string append_suffix_to_filename(const std::string & filename,
		const std::string & sfx) {

	size_t pos = filename.find_last_of(".");
	if (pos == std::string::npos) {
		return filename + sfx;
	}
	return filename.substr(0, pos) + sfx
			+ filename.substr(pos, filename.length());
}


void append_suffix_to_filename(char* out, const std::string & filename,
		const std::string & sfx) {
	strcpy(out, append_suffix_to_filename(filename, sfx).c_str());
}



