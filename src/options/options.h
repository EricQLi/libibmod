/*
 * options.h
 *
 *  Created on: 15 Feb 2013
 *      Author: s12kb2
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <map>
#include <vector>
#include <string>

//#include "anyoption/anyoption.h"
#include "VariableType.h"
#include "string_ext.h"

//{
#define OptionsMapTypes std::string, VariableTypeOption*

class Options {
private:
	unsigned int n_options;
	std::map<OptionsMapTypes> option_list;
	void _setup_option(VariableTypeOption* opt, unsigned int flags);
	void _update_option_list(void);

	bool help_requested;
	std::string help_option_name;

	// environment symbol map
	std::map<std::string, std::string> envSymbols;

	void symbolExpand(std::map<std::string, std::string>& symbols, std::string & s);
	void envSymbolExpand(std::string& s);
//	char file_delimiter_char;


public:
	static VariableTypeOption empty_option;
    bool use_locale;

	std::vector<std::string> get_option_names(void);

	bool processFile(const char *filename, char file_delimiter_char = ':', char commentchar = '#');
	void processCommandArgs(int _argc, char **_argv);

	void setHelpOptionName(const std::string & s);
	bool helpRequested(void);

	void printOptions(void);

	bool saveFile(const char *filename, char file_delimiter_char = ':', char commentchar = '#');
	unsigned int get_option_count(void) const;

	bool is_option_defined(const std::string & name);

	template<class T>
	VariableTypeOption* add_option(T val,
		const std::string & name,
		const std::string & descr,
		unsigned int flags = VariableTypeOption::OPT_DEFAULT) {
		VariableTypeOption* opt = new VariableTypeOption(val, name, descr);
		_setup_option(opt, flags);
		return opt;
	}

	Options(void);
	~Options();

	VariableTypeOption* get_option(std::string name);

	std::string to_string(void);
};

#endif /* OPTIONS_H_ */
