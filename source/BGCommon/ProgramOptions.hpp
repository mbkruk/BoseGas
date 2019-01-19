#ifndef PROGRAM_OPTIONS_HPP_
#define PROGRAM_OPTIONS_HPP_

#include <cinttypes>
#include <iostream>
#include <functional>
#include <memory>
#include <vector>

#include "Trie.hpp"

class ProgramOptions
{
public:

	typedef std::function<int(const char *[])> Callback;

	struct Option
	{
		std::string shortName, longName, description;
		int arguments;
		Callback callback;
	};

private:

	void printDescriptions(Option &opt, std::ostream &os);

	Trie<Option> trie;

	std::vector<std::shared_ptr<Option> > options;

public:

	std::shared_ptr<Option> addOption(const char *shortName, const char *longName, const char *description, int arguments, Callback callback);

	int parseOptions(int argc, const char *argv[]);

	void printDescriptions(std::ostream &os);

	ProgramOptions();
	~ProgramOptions();
};

#endif
