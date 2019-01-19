#ifndef CONFIG_FILE_HPP_
#define CONFIG_FILE_HPP_

#include <fstream>
#include <memory>
#include <functional>

#include "Trie.hpp"

class ConfigFile
{
public:

	typedef std::function<int32_t(std::istream&, const char *,void*)> Callback;

	struct Attribute
	{
		std::string name;
		Callback callback;
	};

private:

	Trie<Attribute> trie;

public:

	void addAttribute(const char *name, Callback callback);

	int32_t parse(const char *fname, void *pData);
};

#endif
