/*
 * options.cpp
 *
 *  Created on: 15 Feb 2013
 *      Author: s12kb2
 */

#include <iostream>
#include <fstream>
#include <map>

#include "options.h"
#include "../console-fun.h"
#include "argvparser.h"

//#include "../defines.h"


#define VECTOR_FOR(_it,_vec,_type) \
	for (std::vector<_type>::iterator _it = _vec.begin(); _it != _vec.end(); ++_it)
#define PMAP_FOR(it, X, ...) for (std::map<__VA_ARGS__>::iterator it = X->begin(); it != X->end(); ++it)
#define MAP_FOR(it, X, ...) for (std::map<__VA_ARGS__>::iterator it = X.begin(); it != X.end(); ++it)


using namespace CommandLineProcessing;


VariableTypeOption Options::empty_option = VariableTypeOption();

Options::Options(void) :
	n_options(0),
	help_requested(false),
	help_option_name("help"),
	use_locale(true) {
}


void Options::setHelpOptionName(const std::string & s) {
	help_option_name = s;
}


// TODO: allow for $NAME.. and $(NAME) or ${NAME} syntax
void Options::symbolExpand(std::map<std::string, std::string>& symbols, std::string & s) {
	bool expanded;
	do {
		expanded = false;
		for (std::map<std::string, std::string>::iterator it = symbols.begin(); it != symbols.end(); ++it) {
			std::string search = "%" + it->first + "%";
			std::string replace = it->second;
			size_t pos = s.find(search);
			if (pos != std::string::npos) {
				expanded = true;
				s.replace(pos, search.length(), replace);
			}
		}
	} while (expanded);
}

void Options::envSymbolExpand(std::string& s) {
	symbolExpand(envSymbols, s);
}


bool Options::is_option_defined(const std::string & name) {
	std::map<OptionsMapTypes>::iterator it = option_list.find(name);
	return it != option_list.end();
}

bool Options::processFile(const char *filename, char separator, char commentchar) {
	FILE* cfgFile = fopen(filename, "r");
	if (!cfgFile) return false;

	char buff[1024];
	while (fgets(buff, 1024, cfgFile)) {
		std::string line = buff;
		if ((line.length() > 2) && (line[0] != commentchar)) {
			std::string name, value;
			strsplit(line, name, value, separator);

			std::map<OptionsMapTypes>::iterator it = option_list.find(name);
			if(it != option_list.end()) {
				envSymbolExpand(value);
				it->second->set_from_str(value);
			} else {
				printf("Unknown option '%s' in config file \n", name.c_str());
			}
		}
	}
	fclose(cfgFile);

//	Config config(filename, 0, _file_delimiter_char);
//
//
//	MAP_FOR(it, option_list, OptionsMapTypes) {
//		std::string name = it->first;
//		if(name == help_option_name) continue;
//		if (config.hasEntry(name))
//			it->second->set_from_str(config.pString(name));
//	}

	return true;
}

/*
bool parseConfigFile(const char* configFile,
		std::map<std::string, std::string> & symbols,
		char** envp, char separator = '=', char commentchar = '#') {
	if(envp != NULL) {
		while (*envp) {
			std::string envEntry = *envp;
			size_t pos = envEntry.find('=');
			if (pos != std::string::npos) {
				std::string name = envEntry.substr(0, pos);
				std::string value = envEntry.substr(pos + 1, std::string::npos);
				//envSymbols[name] = value;
			}
			++envp;
		}
	}
	FILE* in = fopen(configFile, "r");
	if (!in) return false;

	char buff[1024];
	while (fgets(buff, 1024, in)) {

		std::string line = buff;
		if ((line.length() > 2) && (line[0] != commentchar)) {
			std::string name;
			std::string value;
			strsplit(line, name, value, separator);
			envSymbolExpand(value);
			std::map<std::string, std::string>::iterator it;
			it = symbols.find(name);
			if (it != symbols.end()) symbols.erase(it);
			symbols.insert(symbols.end(), std::pair<std::string, std::string>(name, value));

		}
	}
	fclose(in);
	return true;
}
*/

void Options::processCommandArgs(int _argc, char **_argv) {
	ArgvParser parser;

	std::map<OptionsMapTypes>::iterator it = option_list.find(help_option_name);

	if(it != option_list.end())
		parser.setHelpOption(std::string(1, it->second->get_short_name()),
				it->second->get_name(),
				it->second->description);


	VariableTypeOption* x;
	ArgvParser::OptionAttributes optattr;
	MAP_FOR(it, option_list, OptionsMapTypes) {
		if(it->first == help_option_name) continue;
		x = it->second;
		optattr = (x->get_type() != VariableType::VT_BOOLEAN)? ArgvParser::OptionRequiresValue :
				ArgvParser::NoOptionAttribute;

		parser.defineOption(x->get_name(), x->description, optattr);
		if(x->has_short_name()) {
			std::string _short_name(1, x->get_short_name());
			parser.defineOptionAlternative(x->get_name(), _short_name);
		}
	}

	ArgvParser::ParserResults results = parser.parse(_argc, _argv);



	switch (results) {
	case ArgvParser::NoParserError:
		break;
	case ArgvParser::ParserHelpRequested:
		help_requested = true;
		return;
		break;
	default:
		printf("Error parsing command line options: \n %s \n", parser.parseErrorDescription(results).c_str());
		exit(334);
		break;
	}

	MAP_FOR(it, option_list, OptionsMapTypes) {
		x = it->second;
		std::string name = x->get_name();
		if (parser.foundOption(name)) {
			if(x->get_type() == VariableType::VT_BOOLEAN && parser.optionValue(name).empty()) {
				x->set(true);
			} else {
				x->set_from_str(parser.optionValue(name));
			}
		}
	}
}

bool Options::helpRequested(void) {
	return help_requested;
}

void Options::printOptions(void) {
	int shnoff = 1;
	int shn_sfx_len = 9;

	std::vector<std::string> opt_names = get_option_names();
	int max_opt_name_len = 0;
	VECTOR_FOR(it, opt_names, std::string) {
		int len = (int) (*it).length() + (get_option(*it)->has_short_name() ? shn_sfx_len : 0);
		max_opt_name_len = std::max(max_opt_name_len, len);
	}
	max_opt_name_len += shnoff;

	Factor type_labels(0, "undef;bool;int;real;string;range;factor", ";");

	VariableTypeOption* x;
	VECTOR_FOR(it, opt_names, std::string) {
		x = get_option(*it);
//		int  opt_name_len = (*it).length();
		bool has_shortname = x->has_short_name();

		con_set_color(7);
		if(has_shortname) {
			printf("%s (or '%c'):", x->get_name().c_str(), x->get_short_name());
		} else {
			printf("%s:", x->get_name().c_str());
		}
		con_set_color(0);

		std::string space(std::max((size_t) 0, max_opt_name_len - x->get_name().length() - (has_shortname ? 9 : 0)), ' ');
		printf("%s %s ", space.c_str(), x->description.c_str());
		if(x->get_type() == VariableType::VT_FACTOR) {
			printf("(");
			Factor ff = x->get_factor();
			int nlev = ff.getNLevels();
			for (int i = 0; i < nlev; ++i) {
				char fflab[128];
				ff.getLevel(i, fflab);
				printf("%s%s", fflab, (i == nlev - 1) ? ") " : ", ");
			}
		}
		printf("(");
		con_set_color(4);
		std::string default_value, current_value;
		default_value = x->get_default_value()->to_string();
		current_value = x->to_string();
		if(default_value == current_value) 	printf("%s", current_value.c_str());
		else printf("%s, %s", default_value.c_str(), current_value.c_str());
		con_set_color(0);
		printf(") [");
		con_set_color(5);
		printf("%s", type_labels.getLevel((size_t) x->get_type()).c_str());
		con_set_color(0);
		printf("] \n");
	}
}


bool Options::saveFile(const char *filename, char file_delimiter_char, char commentchar) {
    std::ofstream outfile;
	outfile.open(filename);
	if (outfile.is_open()) {
		MAP_FOR(it, option_list, OptionsMapTypes) {
			VariableTypeOption* x = it->second;
				outfile << "# " << x->description << " (" << x->get_default_value()->to_string(false) << ")" << std::endl <<
			x->get_name() << ": " << x->to_string(false) << std::endl;
		}
		outfile.close();
		return true;
	}
	return false;
}

#define FLAG_HAS(f, X)  ((f & (X)) == (X))


void Options::_setup_option(VariableTypeOption* opt, unsigned int flags) {
	std::map<OptionsMapTypes>::iterator it;
	it = option_list.find(opt->get_name());
	if (it != option_list.end()) {
		option_list.erase(it);
		std::cerr << "Replacing element '" << opt->get_name() << "'"
				<< std::endl;
	} else {
		n_options++;
	}
	option_list.insert(option_list.end(), std::pair<OptionsMapTypes>(opt->get_name(), opt));
}

Options::~Options() {
	//std::cout << "Destroying '" << name << "' with value " << to_string() << std::endl;
	MAP_FOR(it1, option_list, OptionsMapTypes) delete it1->second;
	option_list.clear();
}

unsigned int Options::get_option_count(void) const {
	return n_options;
}

std::string Options::to_string(void) {
	return "<Not implemented yet>";
}

std::vector<std::string> Options::get_option_names(void) {
	std::vector<std::string> res;
	MAP_FOR(it, option_list, OptionsMapTypes) {
		res.push_back(it->first);
	}
	return res;
}

VariableTypeOption* Options::get_option(std::string name) {
	std::map<OptionsMapTypes>::iterator it;
	it = option_list.find(name);
	if (it != option_list.end()) {
		return it->second;
	}
	return &empty_option;
}

