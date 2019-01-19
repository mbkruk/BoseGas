#include "ProgramOptions.hpp"
#include <iostream>

std::shared_ptr<ProgramOptions::Option> ProgramOptions::addOption(const char *shortName, const char *longName, const char *description, int arguments, Callback callback)
{
	auto pOpt = std::make_shared<Option>();
	if (description)
		pOpt->description = description;
	pOpt->arguments = arguments;
	pOpt->callback = callback;

	if (shortName)
	{
		pOpt->shortName = shortName;
		trie.set(shortName)->pPayload = pOpt;
	}

	if (longName)
	{
		pOpt->longName = longName;
		trie.set(longName)->pPayload = pOpt;
	}

	options.push_back(pOpt);

	return pOpt;
}

int ProgramOptions::parseOptions(int argc, const char *argv[])
{
	int r;
	for (int i=1;i<argc; )
	{
		auto p = trie.find(argv[i]);
		if (p->pPayload)
		{
			if (argc-i<p->pPayload->arguments)
			{
				std::cerr << "option `" << argv[i] << "` has to have " << p->pPayload->arguments-1 << " arguments" << std::endl;
				return 2;
			}
			if (r=p->pPayload->callback(argv+i))
				return r;
			i += p->pPayload->arguments;
		}
		else
			++i;
	}
	return 0;
}

void ProgramOptions::printDescriptions(Option &opt, std::ostream &os)
{
	if (opt.shortName.empty() && opt.longName.empty())
		return;
	size_t s = opt.shortName.length()+opt.longName.length();
	if (!opt.shortName.empty() && !opt.longName.empty())
		s += 2;
	while (s++<32)
		os << ' ';
	if (!opt.shortName.empty())
	{
		os << opt.shortName;
		if (!opt.longName.empty())
			os << ", ";
	}
	if (!opt.longName.empty())
		os << opt.longName;
	os << " -- " << opt.description << std::endl;
}

void ProgramOptions::printDescriptions(std::ostream &os)
{
	for (auto &pOpt: options)
	if(pOpt)
		printDescriptions(*pOpt,os);
}

ProgramOptions::ProgramOptions()
{
}

ProgramOptions::~ProgramOptions()
{
}
