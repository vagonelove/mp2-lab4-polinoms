#include "monom.h"

MonomList::MonomList()
{
	a = new Monom;
	a->next = a;
}

MonomList::MonomList(const MonomList& poly)
{
	a = new Monom(poly.a->cf, poly.a->deg);
	Monom* p = a;
	Monom* q = poly.a->next;
	while (q != poly.a)
	{
		p->next = new Monom(q->cf, q->deg);
		p = p->next;
		q = q->next;
	}
	p->next = a;
}

MonomList& MonomList::operator=(const MonomList& poly)
{
	if (this != &poly)
	{
		a->cf = poly.a->cf;
		a->deg = poly.a->deg;
		Monom* p = a->next;
		Monom* q = a->next;
		while (p != a)
		{
			q = p->next;
			delete p;
			p = q;
		}
		q = poly.a->next;
		while (q != poly.a)
		{
			p->next = new Monom(q->cf, q->deg);
			p = p->next;
			q = q->next;
		}
		p->next = a;
	}
	return *this;
}

MonomList::~MonomList()
{
	Monom* p = a->next;
	Monom* q = a->next;
	while (p != a)
	{
		q = p->next;
		delete p;
		p = q;
	}
	delete a;
}

void MonomList::Insert(const Monom& monom)
{
	if (abs(monom.cf) > 1e-10)
	{
		Monom* p = a;
		while (monom.deg < p->next->deg)
		{
			p = p->next;
		}
		if (monom.deg == p->next->deg)
		{
			p->next->cf += monom.cf;
			if (p->next->cf < 1e-10 && p->next != a)
			{
				Monom* tmp = p->next;
				p->next = tmp->next;
				delete tmp;
			}
		}
		else
		{
			Monom* tmp = p->next;
			p->next = new Monom(monom.cf, monom.deg, tmp);
		}
	}
}

MonomList::Iterator MonomList::begin()
{
	return MonomList::Iterator(a);
}

MonomList::Iterator::Iterator(Monom* ptr_)
{
	this->ptr = ptr_;
}

MonomList::Iterator& MonomList::Iterator::operator++()
{
	ptr = ptr->next;
	return *this;
}

Monom& MonomList::Iterator::operator*()
{
	return *ptr;
}