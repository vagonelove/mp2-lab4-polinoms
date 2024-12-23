#include "polinom.h"

bool check_err_add_deg(int deg1, int deg2)
{
	int deg1_z = deg1 % 10, deg2_z = deg2 % 10;
	if (deg1_z + deg2_z > 9)
		return true;
	int deg1_y = deg1 % 100 / 10, deg2_y = deg2 % 100 / 10;
	if (deg1_y + deg2_y > 9)
		return true;
	int deg1_x = deg1 / 100, deg2_x = deg2 / 100;
	if (deg1_x + deg2_x > 9)
		return true;
	return false;
}

Polinom& Polinom::operator+=(const Monom& monom)
{
	Insert(monom);
	return *this;
}

Polinom& Polinom::operator-=(const Monom& monom)
{
	Insert(Monom(monom.cf * -1, monom.deg));
	return *this;
}

Polinom Polinom::operator+(const Monom& monom) const
{
	Polinom res(*this);
	res.Insert(monom);
	return res;
}

Polinom Polinom::operator-(const Monom& monom) const
{
	Polinom res(*this);
	res.Insert(Monom(monom.cf * -1, monom.deg));
	return res;
}

Polinom Polinom::operator*(const Monom& monom) const
{
	Polinom res;
	if (monom.deg == 0 || abs(monom.cf) < 1e-10)
		res *= (monom.cf);
	else
	{
		Monom* pThis = a->next;
		Monom* pRes = res.a;

		while (pThis != a)
		{
			if (check_err_add_deg(pThis->deg, monom.deg))
				throw std::exception("Wrong degree!");
			pRes->next = new Monom(pThis->cf * monom.cf, pThis->deg + monom.deg, res.a);
			pThis = pThis->next;
			pRes = pRes->next;
		}
		if (abs(a->cf) > 1e-10)
			pRes->next = new Monom(a->cf * monom.cf, a->deg + monom.deg, res.a);
	}
	return res;
}

Polinom& Polinom::operator*=(const Monom& monom)
{
	if (monom.deg == 0 || abs(monom.cf) < 1e-10)
		operator*=(monom.cf);
	else
	{
		Monom* pThis = a->next;
		while (pThis != a)
		{
			if (check_err_add_deg(pThis->deg, monom.deg))
				throw std::exception("Wrong degree!");
			pThis->cf *= monom.cf;
			pThis->deg += monom.deg;
			pThis = pThis->next;
		}
		if (abs(a->cf) > 1e-10)
		{
			a->cf *= monom.cf;
			a->deg += monom.deg;
			Monom* tmp = a->next;
			a->next = new Monom(0.0, 0, tmp);
			a = a->next;
		}
	}
	return *this;
}

Polinom& Polinom::operator+=(double scalar)
{
	a->cf += scalar;
	return *this;
}

Polinom& Polinom::operator-=(double scalar)
{
	a->cf -= scalar;
	return *this;
}

Polinom Polinom::operator+(double scalar) const
{
	Polinom res(*this);
	res.a->cf += scalar;
	return res;
}

Polinom Polinom::operator-(double scalar) const
{
	Polinom res(*this);
	res.a->cf -= scalar;
	return res;
}

Polinom Polinom::operator*(double scalar) const
{
	Polinom res;
	if (abs(scalar) > 1e-10)
	{
		res.a->cf *= scalar;
		Monom* pThis = a->next;
		Monom* pRes = res.a;
		while (pThis != a)
		{
			pRes->next = new Monom(pThis->cf * scalar, pThis->deg);
			pThis = pThis->next;
			pRes = pRes->next;
		}
	}
	return res;
}

Polinom& Polinom::operator*=(double scalar)
{
	a->cf *= scalar;
	if (abs(scalar) < 1e-10)
	{
		Monom* p = a->next;
		Monom* q = a->next;
		while (p != a)
		{
			q = p->next;
			delete p;
			p = q;
		}
		a->next = a;
	}
	else
	{
		Monom* pThis = a->next;
		while (pThis != a)
		{
			pThis->cf *= scalar;
			pThis = pThis->next;
		}
	}
	return *this;
}

Polinom& Polinom::operator+=(const Polinom& poly)
{
	a->cf += poly.a->cf;
	Monom* pThis = a;
	Monom* pPoly = poly.a->next;
	while (pThis->next != a || pPoly != poly.a)
	{
		if (pPoly->deg > pThis->next->deg)
		{
			Monom* tmp = pThis->next;
			pThis->next = new Monom(pPoly->cf, pPoly->deg, tmp);
			pPoly = pPoly->next;
			pThis = pThis->next;
		}
		else if (pPoly->deg < pThis->next->deg)
		{
			pThis = pThis->next;
		}
		else
		{
			pThis->next->cf += pPoly->cf;
			if (abs(pThis->next->cf) < 1e-10)
			{
				Monom* tmp = pThis->next;
				pThis->next = pThis->next->next;
				delete tmp;
			}
			else
				pThis = pThis->next;

			pPoly = pPoly->next;
		}
	}
	while (pPoly != poly.a)
	{
		pThis->next = new Monom(pPoly->cf, pPoly->deg, a);
		pPoly = pPoly->next;
		pThis = pThis->next;
	}
	return *this;
}

Polinom& Polinom::operator-=(const Polinom& poly)
{
	a->cf -= poly.a->cf;
	if (this == &poly)
	{
		Monom* p = a->next;
		Monom* q = a->next;
		while (p != a)
		{
			q = p->next;
			delete p;
			p = q;
		}
		a->next = a;
	}
	else
	{
		Monom* pThis = a;
		Monom* pPoly = poly.a->next;
		while (pThis->next != a || pPoly != poly.a)
		{
			if (pPoly->deg > pThis->next->deg)
			{
				Monom* tmp = pThis->next;
				pThis->next = new Monom(-pPoly->cf, pPoly->deg, tmp);
				pPoly = pPoly->next;
				pThis = pThis->next;
			}
			else if (pPoly->deg < pThis->next->deg)
			{
				pThis = pThis->next;
			}
			else
			{
				pThis->next->cf -= pPoly->cf;
				if (abs(pThis->next->cf) < 1e-10)
				{
					Monom* tmp = pThis->next;
					pThis->next = pThis->next->next;
					delete tmp;
				}
				else
					pThis = pThis->next;

				pPoly = pPoly->next;
			}
		}
		while (pPoly != poly.a)
		{
			pThis->next = new Monom(-pPoly->cf, pPoly->deg, a);
			pPoly = pPoly->next;
			pThis = pThis->next;
		}
	}
	return *this;
}

Polinom Polinom::operator+(const Polinom& poly) const
{
	Polinom res;
	res.a->cf = a->cf + poly.a->cf;
	Monom* pThis = a->next;
	Monom* pPoly = poly.a->next;
	Monom* pRes = res.a;
	while (pThis != a || pPoly != poly.a)
	{
		if (pThis->deg > pPoly->deg)
		{
			pRes->next = new Monom(pThis->cf, pThis->deg);
			pRes = pRes->next;
			pThis = pThis->next;
		}
		else if (pThis->deg < pPoly->deg)
		{
			pRes->next = new Monom(pPoly->cf, pPoly->deg);
			pRes = pRes->next;
			pPoly = pPoly->next;
		}
		else
		{
			if (abs(pPoly->cf + pThis->cf) > 1e-10)
			{
				pRes->next = new Monom(pPoly->cf + pThis->cf, pPoly->deg);
				pRes = pRes->next;
			}
			pThis = pThis->next;
			pPoly = pPoly->next;
		}
	}
	while (pThis != a)
	{
		pRes->next = new Monom(pThis->cf, pThis->deg);
		pRes = pRes->next;
		pThis = pThis->next;
	}
	while (pPoly != poly.a)
	{
		pRes->next = new Monom(pPoly->cf, pPoly->deg);
		pRes = pRes->next;
		pPoly = pPoly->next;
	}
	pRes->next = res.a;
	return res;
}

Polinom Polinom::operator-(const Polinom& poly) const
{
	Polinom res;
	res.a->cf = a->cf - poly.a->cf;
	Monom* pThis = a->next;
	Monom* pPoly = poly.a->next;
	Monom* pRes = res.a;
	while (pThis != a || pPoly != poly.a)
	{
		if (pThis->deg > pPoly->deg)
		{
			pRes->next = new Monom(pThis->cf, pThis->deg);
			pRes = pRes->next;
			pThis = pThis->next;
		}
		else if (pThis->deg < pPoly->deg)
		{
			pRes->next = new Monom(-pPoly->cf, pPoly->deg);
			pRes = pRes->next;
			pPoly = pPoly->next;
		}
		else
		{
			if (abs(pThis->cf - pPoly->cf) > 1e-10)
			{
				pRes->next = new Monom(pThis->cf - pPoly->cf, pPoly->deg);
				pRes = pRes->next;
			}
			pThis = pThis->next;
			pPoly = pPoly->next;
		}
	}
	while (pThis != a)
	{
		pRes->next = new Monom(pThis->cf, pThis->deg);
		pRes = pRes->next;
		pThis = pThis->next;
	}
	while (pPoly != poly.a)
	{
		pRes->next = new Monom(-pPoly->cf, pPoly->deg);
		pRes = pRes->next;
		pPoly = pPoly->next;
	}
	pRes->next = res.a;
	return res;
}

Polinom Polinom::operator*(const Polinom& poly) const
{
	Polinom res;
	Monom* pThis = a;
	Monom* pPoly = poly.a;
	do
	{
		res += poly * (*pThis);
		pThis = pThis->next;
	} while (pThis != a);
	return res;
}

bool Polinom::operator==(const Polinom& poly) const
{
	Monom* pThis = a;
	Monom* pPoly = poly.a;
	do
	{
		if (*pThis != *pPoly)
			return false;
		pThis = pThis->next;
		pPoly = pPoly->next;
	} while (pThis != a || pPoly != poly.a);
	if (pThis != a || pPoly != poly.a)
		return false;

	return true;
}

void Polinom::str_to_poly(const std::string& input)
{
	parser(input);
	if (isCorrect())
		converter();
	else
		throw std::exception("Error!");
}

void Polinom::parser(const std::string& input)
{
	Lexs.clear();
	int state = 0;
	std::string lexema;
	for (size_t i = 0; i < input.length(); i++)
	{
		if (input[i] != ' ')
		{
			char current = input[i];

			switch (state)
			{
			case 0: {
				if (current == '+' || current == '-')
					state = 1;
				else if (std::isdigit(current) || current == '.')
					state = 2;
				else if (current == 'x' || current == 'y' || current == 'z')
					state = 3;
				else
					throw std::exception("Wrong symbol!!!");
				lexema.append(1, current);
				break;
			}
			case 1: {
				if (std::isdigit(current) || current == '.')
					state = 2;
				else if (current == 'x' || current == 'y' || current == 'z')
					state = 3;
				else
					throw std::exception("Wrong symbol!!!");

				lexema.append(1, current);
				break;
			}
			case 2: {
				if (current == '+' || current == '-')
				{
					state = 1;
					Lexs.push_back(lexema);
					lexema.clear();
				}
				else if (std::isdigit(current) || current == '.')
					state = 2;
				else if (current == 'x' || current == 'y' || current == 'z')
					state = 3;
				else
					throw std::exception("Wrong symbol!!!");

				lexema.append(1, current);
				break;
			}
			case 3: {
				if (current == '+' || current == '-')
				{
					state = 1;
					Lexs.push_back(lexema);
					lexema.clear();
				}
				else if (current == 'x' || current == 'y' || current == 'z')
					state = 3;
				else if (current == '^')
					state = 4;
				else
					throw std::exception("Wrong symbol!!!");

				lexema.append(1, current);
				break;
			}
			case 4: {
				if (std::isdigit(current) || current == '.')
				{
					state = 5;
					lexema.append(1, current);
				}
				else
					throw std::exception("Wrong symbol!!!");

				break;
			}
			case 5: {
				if (current == '+' || current == '-')
				{
					state = 1;
					Lexs.push_back(lexema);
					lexema.clear();
				}
				else if (current == 'x' || current == 'y' || current == 'z')
					state = 3;
				else
					throw std::exception("Wrong symbol!!!");

				lexema.append(1, current);
				break;
			}
			default:
				break;
			}
		}
	}
	if (state == 2 || state == 3 || state == 5)
		Lexs.push_back(lexema);
	else
		throw std::exception("Wrong symbol!!!");
}

bool Polinom::check_vars() const
{
	for (size_t i = 0; i < Lexs.size(); i++)
	{
		size_t pX = Lexs[i].find_first_of('x');
		pX = Lexs[i].find_first_of('x', pX + 1);
		size_t pY = Lexs[i].find_first_of('y');
		pY = Lexs[i].find_first_of('y', pY + 1);
		size_t pZ = Lexs[i].find_first_of('z');
		pZ = Lexs[i].find_first_of('z', pZ + 1);
		if (pX != Lexs[i].npos || pY != Lexs[i].npos || pZ != Lexs[i].npos)
			return false;
	}
	return true;
}

bool Polinom::check_points() const
{
	for (size_t i = 0; i < Lexs.size(); i++)
	{
		size_t p = Lexs[i].find_first_of('.');
		p = Lexs[i].find_first_of('.', p + 1);
		if (p != Lexs[i].npos)
			return false;
	}
	return true;
}

bool Polinom::isCorrect() const
{
	return check_points() && check_vars();
}

void Polinom::converter()
{
	for (size_t i = 0; i < Lexs.size(); i++)
	{
		int state = 0;
		std::string cf;
		int degree = 0;
		for (size_t j = 0; j < Lexs[i].length(); j++)
		{
			const char current = Lexs[i][j];

			switch (state)
			{
			case 0: {
				if (std::isdigit(current) || current == '.')
				{
					state = 2;
					cf.append(1, current);
				}
				else if (current == '-' || current == '+')
				{
					state = 1;
					cf.append(1, current);
				}
				else
				{
					state = 3;
					cf.append("1");
				}

				break;
			}
			case 1: {
				if (std::isdigit(current) || current == '.')
				{
					state = 2;
					cf.append(1, current);
				}
				else
				{
					state = 3;
					cf.append("1");
				}
				break;
			}
			case 2: {
				if (std::isdigit(current) || current == '.')
				{
					state = 2;
					cf.append(1, current);
				}
				else
				{
					state = 3;
				}
				break;
			}
			case 3: {
				if (current == '^')
				{
					if (Lexs[i][j - 1] == 'x')
						degree += std::stoi(&Lexs[i][j + 1]) * 100;
					else if (Lexs[i][j - 1] == 'y')
						degree += std::stoi(&Lexs[i][j + 1]) * 10;
					else
						degree += std::stoi(&Lexs[i][j + 1]);
					j++;
				}
				else
				{
					if (Lexs[i][j - 1] == 'x')
						degree += 100;
					else if (Lexs[i][j - 1] == 'y')
						degree += 10;
					else if (Lexs[i][j - 1] == 'z')
						degree += 1;
				}
				break;
			}
			default:
				break;
			}
		}
		if (Lexs[i].back() == 'x')
			degree += 100;
		else if (Lexs[i].back() == 'y')
			degree += 10;
		else if (Lexs[i].back() == 'z')
			degree += 1;

		Insert(Monom(std::stod(cf), degree));
	}
}

void Polinom::Print() const
{
	Monom* p = a->next;
	while (p != a)
	{
		std::cout << p->cf << "x^" << p->deg / 100 << "y^" << (p->deg % 100) / 10 << "z^" << p->deg % 10 << " + ";
		p = p->next;
	}
	std::cout << a->cf;
	std::cout << std::endl;
}