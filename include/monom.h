// класс или структура Монома (double коэффициент, целая свернутая степень)

#pragma once
#include <iostream>

struct Monom
{
	int deg;
	double cf;
	Monom* next;
	Monom(double cf_ = 0.0, int deg_ = 0, Monom* next_ = nullptr)
	{
		if (deg_ < 0 || deg_ > 999)
			throw std::exception("Wrong monom's degree!");
		deg = deg_;
		cf = cf_;
		next = next_;
	}
	Monom(const Monom& monom2)
	{
		deg = monom2.deg;
		cf = monom2.cf;
		next = monom2.next;
	}
	Monom& operator=(double scalar)
	{
		deg = 0;
		cf = scalar;
		next = nullptr;
	}
	bool operator==(const Monom& m) const
	{
		return (abs(cf) - abs(m.cf) < 1e-10) && (deg == m.deg);
	}
	bool operator!=(const Monom& m) const
	{
		return (abs(cf) - abs(m.cf) > 1e-10) || (deg != m.deg);
	}
	friend std::ostream& operator<<(std::ostream& os, const Monom& m)
	{
		os << m.cf << "x^" << m.deg / 100 << "y^" << (m.deg % 100) / 10 << "z^" << m.deg % 10 << std::endl;
		return os;
	}
};

class MonomList
{
protected:
	Monom* a;
public:
	MonomList();
	MonomList(const MonomList& poly);
	MonomList& operator=(const MonomList& poly);
	~MonomList();
	void Insert(const Monom& monom);
	class Iterator
	{
		Monom* ptr;
	public:
		Iterator(Monom* _ptr = nullptr);
		Iterator& operator++();
		Monom& operator*();
	};
	Iterator begin();
};