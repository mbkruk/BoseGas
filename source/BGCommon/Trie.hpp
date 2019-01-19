#ifndef TRIE_HPP_
#define TRIE_HPP_

#include <cinttypes>
#include <memory>

template<typename Payload> class Trie
{
public:

	static inline uint32_t getCharIndex(char code)
	{
		if ('A'<=code && code<='Z')
			return code-'A';
		if ('a'<=code && code<='z')
			return code-'a'+26;
		if (code=='-')
			return 52;
		return 0xffffffff;
	}

	struct CharLeaf
	{
		static constexpr uint32_t kCount = 53;
		char code;
		CharLeaf* pNext[kCount];

		std::shared_ptr<Payload> pPayload;

		void initialize()
		{
			code = 0;
			for (uint32_t i=0;i<kCount;++i)
				pNext[i] = nullptr;
		}

		void release()
		{
			for (uint32_t i=0;i<kCount;++i)
			if (pNext[i])
			{
				pNext[i]->release();
				delete pNext[i];
				pNext[i] = nullptr;
			}
			pPayload.reset();
		}

	};

protected:

	CharLeaf trieRoot;

public:

	inline std::shared_ptr<Payload> root()
	{
		return trieRoot.pPayload;
	}

	CharLeaf* set(const char *word)
	{
		CharLeaf *p = &trieRoot;
		const char *c = word;
		while(*c)
		{
			if (*c==' ')
				break;
			uint32_t i = getCharIndex(*c);
			if (!p->pNext[i])
			{
				p->pNext[i] = new CharLeaf;
				p->pNext[i]->initialize();
			}
			p = p->pNext[i];
			++c;
		}
		return p;
	}

	CharLeaf* find(const char *word)
	{
		CharLeaf *p = &trieRoot;
		const char *c = word;
		while (*c)
		{
			uint32_t j = getCharIndex(*c);
			if (j==0xffffffff || !p->pNext[j])
			{
				p = &trieRoot;
				break;
			}
			p = p->pNext[j];
			++c;
		}
		return p;
	}

	Trie()
	{
		trieRoot.initialize();
	}

	~Trie()
	{
		trieRoot.release();
	}
};

#endif
