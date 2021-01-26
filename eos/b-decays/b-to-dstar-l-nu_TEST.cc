/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Christoph Bobeth
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
#include <eos/observable.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <iomanip>

using namespace test;
using namespace eos;

class BToVectorLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToVectorLeptonNeutrinoTest() :
            TestCase("b_to_dstar_l_nu_test")
        {
        }

        virtual void run() const
        {
            // comparison with Martin Jung in 3/2/1 model
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["mass::e"].set(0.000001);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "e"       },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    }
                };

                BToVectorLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                // Christoph Bobeth: Adjusted test case because increased number of integration points
                //                   in numerical integration from 256 -> 2048
                TEST_CHECK_NEARLY_EQUAL(/*33.323*/33.3247,       d.integrated_branching_ratio(0.001, 10.689), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.546,       d.integrated_f_L(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.409302220, d.integrated_S1c(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.255523335, d.integrated_S1s(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.409302220, d.integrated_S2c(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.085174445, d.integrated_S2s(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.134468151, d.integrated_S3 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.231808464, d.integrated_S4 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.165381861, d.integrated_S5 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S6c(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.200153929, d.integrated_S6s(0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0,         d.integrated_S7 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S8 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0,         d.integrated_S9 (0.001, 10.689),  eps);
                
                // Some numerics with physical electron mass
                // for comparison with Martin:
                //    rho^2=-xi', cxi=1/2 x'', dxi=1/6 x'''
                //
                Parameters p2 = Parameters::Defaults();

                p2["B(*)->D(*)::xi'(1)@HQET"].set(-1.14);
                p2["B(*)->D(*)::xi''(1)@HQET"].set(2* 0.94);
                p2["B(*)->D(*)::xi'''(1)@HQET"].set(6* (-0.548333));
                p2["B(*)->D(*)::chi_2(1)@HQET"].set(-0.06);
                p2["B(*)->D(*)::chi_2'(1)@HQET"].set(0.0);
                p2["B(*)->D(*)::chi_2''(1)@HQET"].set(0.06);
                p2["B(*)->D(*)::chi_3'(1)@HQET"].set(0.04);
                p2["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.05);
                p2["B(*)->D(*)::eta(1)@HQET"].set(0.6);
                p2["B(*)->D(*)::eta'(1)@HQET"].set(-0.02);
                p2["B(*)->D(*)::eta''(1)@HQET"].set(-0.04);
                p2["B(*)->D(*)::l_1(1)@HQET"].set(0.12);
                p2["B(*)->D(*)::l_1'(1)@HQET"].set(-5.78);
                p2["B(*)->D(*)::l_2(1)@HQET"].set(-1.89);
                p2["B(*)->D(*)::l_2'(1)@HQET"].set(-3.14);
                p2["B(*)->D(*)::l_3(1)@HQET"].set(0.86);
                p2["B(*)->D(*)::l_3'(1)@HQET"].set(0.06);
                p2["B(*)->D(*)::l_4(1)@HQET"].set(-2.02);
                p2["B(*)->D(*)::l_4'(1)@HQET"].set(-0.05);
                p2["B(*)->D(*)::l_5(1)@HQET"].set(3.79);
                p2["B(*)->D(*)::l_5'(1)@HQET"].set(-1.4);
                p2["B(*)->D(*)::l_6(1)@HQET"].set(3.53);
                p2["B(*)->D(*)::l_6'(1)@HQET"].set(0.04);

                p2["CKM::abs(V_cb)"].set(1.0);
                p2["mass::e"].set(0.00000000000001);
                p2["B(*)->D(*)::a@HQET"].set(1.000);
                p2["mass::B_d"].set(5.27942);
                p2["mass::D_u^*"].set(2.01000);
                p2["mass::D_d^*"].set(2.01000);
                p2["life_time::B_d"].set(1.520e-12);
                
                Options o2mu{
                    { "l",             "mu"      },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    }
                };
                
                Options o2e{
                    { "l",             "e"      },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    }
                };
                
                Kinematics k2mu
                {
                    { "q2_min",     p2["mass::mu"]* p2["mass::mu"] }, { "q2_max",    10.689 }, 
                    { "q2_e_min",   p2["mass::e"] * p2["mass::e"]  }, { "q2_e_max",  10.689 },
                    { "q2_mu_min",  p2["mass::mu"]* p2["mass::mu"] }, { "q2_mu_max", 10.689 },
                };

                Kinematics k2e
                {
                    { "q2", 0.0000000000001 },
                    { "q2_min",     p2["mass::e"] * p2["mass::e"]  }, { "q2_max",    10.689 }, 
                    { "q2_e_min",   p2["mass::e"] * p2["mass::e"]  }, { "q2_e_max",  10.689 },
                    { "q2_mu_min",  p2["mass::mu"]* p2["mass::mu"] }, { "q2_mu_max", 10.689 },
                };

                Kinematics kBelle  // w = 1.5 corresponds to q2_min = 0.0775
                {
                    { "q2_min",     0.0775 }, { "q2_max",    10.689 }, 
                    { "q2_e_min",   0.0775 }, { "q2_e_max",  10.689 },
                    { "q2_mu_min",  0.0775 }, { "q2_mu_max", 10.689 },
                };

                auto obs_dBr_e             = Observable::make("B->D^*lnu::dBR/dq2", p2, k2e, o2e);

                auto obs_Br_CPave_mu       = Observable::make("B->D^*lnu::BR^CPave", p2, k2mu, o2mu);
                auto obs_Belle_Br_CPave_mu = Observable::make("B->D^*lnu::BR^CPave", p2, kBelle, o2mu);
                auto obs_Br_CPave_e        = Observable::make("B->D^*lnu::BR^CPave", p2, k2e, o2e);
                auto obs_Belle_Br_CPave_e  = Observable::make("B->D^*lnu::BR^CPave", p2, kBelle, o2e);
 
                auto obs_Afb_CPave_mu       = Observable::make("B->D^*lnu::A_FB^CPave", p2, k2mu, o2mu);
                auto obs_Belle_Afb_CPave_mu = Observable::make("B->D^*lnu::A_FB^CPave", p2, kBelle, o2mu);
                auto obs_Afb_CPave_e        = Observable::make("B->D^*lnu::A_FB^CPave", p2, k2e, o2e);
                auto obs_Belle_Afb_CPave_e  = Observable::make("B->D^*lnu::A_FB^CPave", p2, kBelle, o2e);

                auto obs_FL_e              = Observable::make("B->D^*lnu::F_L",            p2, k2e, o2e);
                auto obs_FL_CPave_e        = Observable::make("B->D^*lnu::F_L^CPave",      p2, k2e, o2e);
                auto obs_FtildeL_e         = Observable::make("B->D^*lnu::Ftilde_L",       p2, k2e, o2e);
                auto obs_FtildeL_CPave_e   = Observable::make("B->D^*lnu::Ftilde_L^CPave", p2, k2e, o2e);
                auto obs_FL_mu             = Observable::make("B->D^*lnu::F_L",            p2, k2mu, o2mu);
                auto obs_FL_CPave_mu       = Observable::make("B->D^*lnu::F_L^CPave",      p2, k2mu, o2mu);
                auto obs_FtildeL_mu        = Observable::make("B->D^*lnu::Ftilde_L",       p2, k2mu, o2mu);
                auto obs_FtildeL_CPave_mu  = Observable::make("B->D^*lnu::Ftilde_L^CPave", p2, k2mu, o2mu);

                auto obs_S3                 = Observable::make("B->D^*lnu::S_3",        p2, k2mu, o2mu);

                // mu-e combined observables
                auto obs_Afb_bar_CPave       = Observable::make("B->D^*lnu::barA_FB^CPave",   p2, k2mu, o2mu);
                auto obs_Afb_del_CPave       = Observable::make("B->D^*lnu::DeltaA_FB^CPave", p2, k2mu, o2mu);
                auto obs_Belle_Afb_bar_CPave = Observable::make("B->D^*lnu::barA_FB^CPave",   p2, kBelle, o2mu);
                auto obs_Belle_Afb_del_CPave = Observable::make("B->D^*lnu::DeltaA_FB^CPave", p2, kBelle, o2mu);

                auto obs_FL_bar_CPave       = Observable::make("B->D^*lnu::barF_L^CPave",    p2, k2mu, o2mu);
                auto obs_FL_del_CPave       = Observable::make("B->D^*lnu::DeltaF_L^CPave",  p2, k2mu, o2mu);
                auto obs_Belle_FL_bar_CPave = Observable::make("B->D^*lnu::barF_L^CPave",    p2, kBelle, o2mu);
                auto obs_Belle_FL_del_CPave = Observable::make("B->D^*lnu::DeltaF_L^CPave",  p2, kBelle, o2mu);

                auto obs_FtildeL_bar_CPave       = Observable::make("B->D^*lnu::barFtilde_L^CPave",   p2, k2mu, o2mu);
                auto obs_FtildeL_del_CPave       = Observable::make("B->D^*lnu::DeltaFtilde_L^CPave", p2, k2mu, o2mu);
                auto obs_Belle_FtildeL_bar_CPave = Observable::make("B->D^*lnu::barFtilde_L^CPave",   p2, kBelle, o2mu);
                auto obs_Belle_FtildeL_del_CPave = Observable::make("B->D^*lnu::DeltaFtilde_L^CPave", p2, kBelle, o2mu);

                auto obs_S3_bar = Observable::make("B->D^*lnu::barS_3",   p2, k2mu, o2mu);
                auto obs_S3_del = Observable::make("B->D^*lnu::DeltaS_3", p2, k2mu, o2mu);
                auto obs_Belle_S3_bar = Observable::make("B->D^*lnu::barS_3",   p2, kBelle, o2mu);
                auto obs_Belle_S3_del = Observable::make("B->D^*lnu::DeltaS_3", p2, kBelle, o2mu);

                std::cout << std::setprecision(20) << std::endl
                << "e   dBr/dq^2         = " << obs_dBr_e->evaluate() << std::endl << std::endl

                << "mu  Br^CPave         = " << obs_Br_CPave_mu->evaluate()       << std::endl
                << "mu  Br^CPave Belle   = " << obs_Belle_Br_CPave_mu->evaluate() << std::endl
                << "e   Br^CPave         = " << obs_Br_CPave_e->evaluate()        << std::endl
                << "e   Br^CPave Belle   = " << obs_Belle_Br_CPave_e->evaluate()  << std::endl

                << "mu  A_FB^CPave       = " << obs_Afb_CPave_mu->evaluate()       << std::endl
                << "mu  A_FB^CPave Belle = " << obs_Belle_Afb_CPave_mu->evaluate() << std::endl
                << "e   A_FB^CPave       = " << obs_Afb_CPave_e->evaluate()        << std::endl
                << "e   A_FB^CPave Belle = " << obs_Belle_Afb_CPave_e->evaluate()  << std::endl
                
                << "e   F_L              = " << obs_FL_e->evaluate()             << std::endl
                << "e   F_L^CPave        = " << obs_FL_CPave_e->evaluate()       << std::endl
                << "e   Ftilde_L         = " << obs_FtildeL_e->evaluate()        << std::endl
                << "e   Ftilde_L^CPave   = " << obs_FtildeL_CPave_e->evaluate()  << std::endl
                << "mu  F_L              = " << obs_FL_mu->evaluate()            << std::endl
                << "mu  F_L^CPave        = " << obs_FL_CPave_mu->evaluate()      << std::endl
                << "mu  Ftilde_L         = " << obs_FtildeL_mu->evaluate()       << std::endl
                << "mu  Ftilde_L^CPave   = " << obs_FtildeL_CPave_mu->evaluate() << std::endl

                << "    S_3              = " << obs_S3->evaluate()                 << std::endl<< std::endl

                << "del A_FB     = " << obs_Afb_del_CPave->evaluate()     << std::endl
                << "del F_L      = " << obs_FL_del_CPave->evaluate()      << std::endl
                << "del Ftilde_L = " << obs_FtildeL_del_CPave->evaluate() << std::endl
                << "del S_3      = " << obs_S3_del->evaluate()            << std::endl
                << "bar A_FB     = " << obs_Afb_bar_CPave->evaluate()     << std::endl
                << "bar F_L      = " << obs_FL_bar_CPave->evaluate()      << std::endl
                << "bar Ftilde_L = " << obs_FtildeL_bar_CPave->evaluate() << std::endl
                << "bar S_3      = " << obs_S3_bar->evaluate()            << std::endl << std::endl
                
                << "bar A_FB Belle     = " << obs_Belle_Afb_bar_CPave->evaluate() << std::endl
                << "del A_FB Belle     = " << obs_Belle_Afb_del_CPave->evaluate() << std::endl
                << "bar F_L Belle      = " << obs_Belle_FL_bar_CPave->evaluate()  << std::endl
                << "del F_L Belle      = " << obs_Belle_FL_del_CPave->evaluate()  << std::endl
                << "bar Ftilde_L Belle = " << obs_Belle_FtildeL_bar_CPave->evaluate() << std::endl
                << "del Ftilde_L Belle = " << obs_Belle_FtildeL_del_CPave->evaluate() << std::endl
                << "bar S_3 Belle      = " << obs_Belle_S3_bar->evaluate() << std::endl
                << "del S_3 Belle      = " << obs_Belle_S3_del->evaluate() << std::endl;

/* Martin:
bar AFB     = -0.196   +- 0.016
del AFB     = -0.00568 +- 0.00031

bar   FL    =  0.546    +- 0.013
delta FL    = -0.000369 +- 0.000045

bar FLtilde =  0.542   +- 0.013
del FLtilde = -0.00713 +- 0.00042

bar S3      = -0.1345   +- 0.0045
del S3      =  0.000311 +- 0.000018

my:

1024 int {1.0466, -44.9628, -1.57842, 9.8461}%

del A_FB     = -0.0055413456334479083143
del F_L      = -0.00025602085231768434426
del Ftilde_L = -0.0067511190567718371014
del S_3      = 0.00034226584179541941211
bar A_FB     = 0.20414019386941184564
bar F_L      = 0.53841012830568191205
bar Ftilde_L = 0.53516257920345466914
bar S_3      = -0.13497880668962561335

2048 int {0.482792, -17.5831, -0.78135, 4.99951}%

del A_FB     = -0.0055099516239523471661
del F_L      = -0.0003156362784531019372
del Ftilde_L = -0.0068045129591636221988
del S_3      = 0.00032480464742329351324
bar A_FB     = 0.20413111557343269209
bar F_L      = 0.53843341137244260519
bar Ftilde_L = 0.53518897303208678995
bar S_3      = -0.13497216630486247313

4096 int {0.222081, -7.87194, -0.394709, 2.50151}%

del A_FB     = -0.0054955545965343910453
del F_L      = -0.00034405146226346161598
del Ftilde_L = -0.006830718515506428723
del S_3      = 0.00031648285315305502152
bar A_FB     = 0.20412605962996527298
bar F_L      = 0.53844587443292823981
bar Ftilde_L = 0.53520254090630592358
bar S_3      = -0.1349685702882043592

5096 int {0.175985, -6.21083, -0.318133, 2.0136}%

del A_FB     = -0.0054930169038115705948
del F_L      = -0.00034943236167350733012
del Ftilde_L = -0.0068359326681268761661
del S_3      = 0.00031490695430208548444
bar A_FB     = 0.20412493715390450655
bar F_L      = 0.53844843038756606823
bar Ftilde_L = 0.53520518023433916177
bar S_3      = -0.13496782883826388688

20384 int {0.0410792, -1.47996, -0.0817706, 0.506238}%

del A_FB     = -0.005485603439460390307
del F_L      = -0.00036572245994570540262
del Ftilde_L = -0.0068520770141193843017
del S_3      = 0.00031013602490290170799
bar A_FB     = 0.20412133286579281499
bar F_L      = 0.53845643605385173913
bar Ftilde_L = 0.53521325877676362293
bar S_3      = -0.13496549689032827368

*/
            }

            // comparison with Martin Jung in 3/2/1 model
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "tau"     },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    }
                };

                BToVectorLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                TEST_CHECK_NEARLY_EQUAL( 8.213,        d.integrated_branching_ratio(3.157, 10.689), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.475,        d.integrated_f_L(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.4325856250, d.integrated_S1c(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2779590234, d.integrated_S1s(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1287773345, d.integrated_S2c(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0471441750, d.integrated_S2s(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0819412032, d.integrated_S3 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1057578408, d.integrated_S4 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2056068494, d.integrated_S5 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2766922602, d.integrated_S6c(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1598442669, d.integrated_S6s(3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0,          d.integrated_S7 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,          d.integrated_S8 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0,          d.integrated_S9 (3.157, 10.689),  eps);
            }

            // SM tests cf. [DSD2014]
            {
                Parameters p1 = Parameters::Defaults();
                /*
                 * for the TEST case below the B->D^* SSE parameters are randomly chosen. However, the correlations together with the EOM conditions among the FFs are respected in this choice. Namely, alpha^A0_0 is correlated with alpha^A12_0, and also alpha^T1_0 should be the same as alpha^T2_0.
                 */
                p1["B->D^*::alpha^A0_0@BSZ2015" ] = +1.0;
                p1["B->D^*::alpha^A0_1@BSZ2015" ] = +0.24;
                p1["B->D^*::alpha^A0_2@BSZ2015" ] = +0.21;
                p1["B->D^*::alpha^A1_0@BSZ2015" ] = +0.5;
                p1["B->D^*::alpha^A1_1@BSZ2015" ] = +0.4;
                p1["B->D^*::alpha^A1_2@BSZ2015" ] = +0.3;
                p1["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p1["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p1["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p1["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p1["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p1["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p1["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p1["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p1["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p1["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p1["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p1["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p1["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p1["mass::B_d"]                   = +5.279;
                p1["mass::D_d^*"]                 = +2.0103;
                // by default, all other couplings are zero in eos
                p1["b->cmunumu::Re{cVL}"]         = +1.0066;  // include Sirlin correction
                p1["b->ctaunutau::Re{cVL}"]       = +1.0066;  // include Sirlin correction

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");
                oo.set("q", "d");

                BToVectorLeptonNeutrino d(p1, oo);

                const double eps = 1e-3;

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 25.4230, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(4.0, 10.68), 0.000494949, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(4.0, 10.68), 0.737489, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(4.0, 10.68), -0.130926, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(4.0, 10.68), 0.00266046, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(4.0, 10.68), 0.230111, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(4.0, 10.68), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(4.0, 10.68), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(4.0, 10.68), 0.0, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p1, k, oo);
                TEST_CHECK_RELATIVE_ERROR(0.379092, obs_RDst->evaluate(), eps);
            }

            // NP tests cf. [DSD2014]
            {
                const double etaEW = 1.0066;
                
                Parameters p3 = Parameters::Defaults();
                /*
                 * for the TEST case below the B->D^* SSE parameters are randomly chosen.
                 * However, the correlations together with the EOM conditions among the FFs
                 * are respected in this choice. Namely, alpha^A0_0 is correlated with alpha^A12_0,
                 * and also alpha^T1_0 should be the same as alpha^T2_0.
                 */
                p3["B->D^*::alpha^A0_0@BSZ2015" ] = +1.0;
                p3["B->D^*::alpha^A0_1@BSZ2015" ] = +0.24;
                p3["B->D^*::alpha^A0_2@BSZ2015" ] = +0.21;
                p3["B->D^*::alpha^A1_0@BSZ2015" ] = +0.5;
                p3["B->D^*::alpha^A1_1@BSZ2015" ] = +0.4;
                p3["B->D^*::alpha^A1_2@BSZ2015" ] = +0.3;
                p3["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p3["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p3["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p3["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p3["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p3["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p3["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p3["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p3["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p3["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p3["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p3["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p3["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p3["mass::B_d"]                   = +5.279;
                p3["mass::D_d^*"]                 = +2.0103;
                // fix scale
                p3["mu"]                          = +4.18;
                // mb(mb)
                p3["mass::b(MSbar)"]              = +4.18;
                // mc(mc)
                p3["mass::c"]                     = +1.275;
                // mu mode
                p3["b->cmunumu::Re{cVL}"]         = +1.0 * etaEW;
                p3["b->cmunumu::Im{cVL}"]         = -2.0 * etaEW;
                p3["b->cmunumu::Re{cVR}"]         = +2.0 * etaEW;
                p3["b->cmunumu::Im{cVR}"]         = -2.0 * etaEW;
                p3["b->cmunumu::Re{cSL}"]         = +3.0 * etaEW;
                p3["b->cmunumu::Im{cSL}"]         = -3.0 * etaEW;
                p3["b->cmunumu::Re{cSR}"]         = +4.0 * etaEW;
                p3["b->cmunumu::Im{cSR}"]         = -4.0 * etaEW;
                p3["b->cmunumu::Re{cT}"]          = +5.0 * etaEW;
                p3["b->cmunumu::Im{cT}"]          = -5.0 * etaEW;
                // tau mode
                p3["b->ctaunutau::Re{cVL}"]       = +1.0 * etaEW;
                p3["b->ctaunutau::Im{cVL}"]       = -5.0 * etaEW;
                p3["b->ctaunutau::Re{cVR}"]       = +2.1 * etaEW;
                p3["b->ctaunutau::Im{cVR}"]       = -6.0 * etaEW;
                p3["b->ctaunutau::Re{cSL}"]       = +3.1 * etaEW;
                p3["b->ctaunutau::Im{cSL}"]       = -7.0 * etaEW;
                p3["b->ctaunutau::Re{cSR}"]       = +4.1 * etaEW;
                p3["b->ctaunutau::Im{cSR}"]       = -8.0 * etaEW;
                p3["b->ctaunutau::Re{cT}"]        = +5.1 * etaEW;
                p3["b->ctaunutau::Im{cT}"]        = -9.0 * etaEW;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");

                BToVectorLeptonNeutrino d(p3, oo);

                const double eps = 1e-3;

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 3431.13, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(4.0, 10.68), 0.0409932, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(4.0, 10.68), 0.50729, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(4.0, 10.68), 0.184031, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(4.0, 10.68), -0.0282197, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(4.0, 10.68), -0.42545, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(4.0, 10.68), 0.0000348895, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(4.0, 10.68), 0.000268975, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(4.0, 10.68), -0.0000320251, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p3, k, oo);
                TEST_CHECK_RELATIVE_ERROR(1.20331, obs_RDst->evaluate(), eps);
            }
        }
} b_to_dstar_l_nu_test;
