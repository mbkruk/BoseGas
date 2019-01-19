#include "ConfigFile.hpp"

void ConfigFile::addAttribute(const char *name, Callback callback)
{
	auto p = std::make_shared<Attribute>();
	p->name = name;
	p->callback = callback;
	trie.set(name)->pPayload = p;
}

int32_t ConfigFile::parse(const char *fname, void *pData)
{
	std::ifstream f(fname);
	if (!f.is_open())
		return 1;
	std::string word;
	int32_t r;
	while (!f.eof())
	{
		if(!(f>>word) || f.eof() || word.empty())
			break;
		auto p = trie.find(word.c_str())->pPayload;
		if (!p)
			return trie.root()->callback(f,word.c_str(),pData);
		if (r=p->callback(f,word.c_str(),pData))
			return r;
	}
	f.close();
	return 0;
}
