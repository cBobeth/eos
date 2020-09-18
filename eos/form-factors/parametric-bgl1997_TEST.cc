/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Christoph Bobeth
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/form-factors/parametric-bgl1997-impl.hh>

#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BGL1997FormFactorsTest :
    public TestCase
{
    public:
        BGL1997FormFactorsTest() :
            TestCase("BGL1997_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            //p["QCD::alpha_s(MZ)"] = 0.1176;

            /* Outer function */
            {
                BGL1997FormFactors<BToDstar> ff(p, Options{ });
                const double m_B(BToDstar::mB);
                const double m_V(BToDstar::mV);
                const double t_0((m_B + m_V)* power_of<2>(std::sqrt(m_B) - std::sqrt(m_V)));

                //TEST_CHECK_NEARLY_EQUAL( 0.0535063, ff._z(1.0, 2.0), eps);

                //_phi(s, t_0, K, a, b, c, chi);
                TEST_CHECK_NEARLY_EQUAL( 0.0331832,  ff._phi(-2.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0324458,  ff._phi(+1.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0316657,  ff._phi(+4.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.488275,   ff._phi(-2.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.470779,   ff._phi(+1.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.452784,   ff._phi(+4.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.00817026, ff._phi(-2.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00822179, ff._phi(+1.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00827232, ff._phi(+4.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
            }

            /* B -> D^* FFs */
            {
                BGL1997FormFactors<BToDstar> ff(p, Options{ });

                // TEST_CHECK_NEARLY_EQUAL( 0.0,     pi.phi(1.0, 2.0), eps);
            }

            /* B -> D FFs*/
            {
                BGL1997FormFactors<BToD> ff(p, Options{ });

                // TEST_CHECK_NEARLY_EQUAL( 0.0045,         pi.f3pi(1.0),        eps);
            }

        }
} BGL1997_form_factor_test;
