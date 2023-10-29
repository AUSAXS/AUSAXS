#include <hist/distance_calculator/DebugManager.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <hist/foxs/FormFactorFoXS.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <constants/Constants.h>
#include <container/Container2D.h>

using namespace hist;
using namespace data::record;

DebugDistanceHistogram::DebugDistanceHistogram(view_ptr<const data::Molecule> protein) : protein(protein) {
    this->q_axis = Axis(settings::axes::qmin, settings::axes::qmax, 100).as_vector();
}

DebugDistanceHistogram::~DebugDistanceHistogram() = default;

ScatteringProfile DebugDistanceHistogram::debye_transform() const {
    static unsigned int counter = 0;
    std::cout << "debye transform " << counter++ << std::endl;

    std::vector<Atom> atoms = protein->get_atoms();
    double Z_exv = protein->get_excluded_volume()*constants::charge::density::water/atoms.size();

    std::vector<double> I;
    I.reserve(q_axis.size());

    container::Container2D<double> aa_distances(atoms.size(), atoms.size());
    container::Container2D<double> aw_distances(atoms.size(), protein->get_waters().size());
    container::Container2D<double> ww_distances(protein->get_waters().size(), protein->get_waters().size());
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = 0; j < atoms.size(); ++j) {
            aa_distances(i, j) = atoms[i].distance(atoms[j]);
        }
        for (unsigned int j = 0; j < protein->get_waters().size(); ++j) {
            aw_distances(i, j) = atoms[i].distance(protein->get_waters()[j]);
        }
    }
    for (unsigned int i = 0; i < protein->get_waters().size(); ++i) {
        for (unsigned int j = 0; j < protein->get_waters().size(); ++j) {
            ww_distances(i, j) = protein->get_waters()[i].distance(protein->get_waters()[j]);
        }
    }

    std::vector<double> count( constants::axes::d_axis.bins, 0);
    std::vector<double> xx_sum(constants::axes::d_axis.bins, 0);
    std::vector<double> aa_sum(constants::axes::d_axis.bins, 0);
    hist::detail::CompactCoordinatesFF cc(protein->get_bodies());
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = i+1; j < atoms.size(); ++j) {
            double fi = form_factor::foxs::storage::atomic::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(i))).evaluate(0);
            double fj = form_factor::foxs::storage::atomic::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(j))).evaluate(0);

            double fix = form_factor::foxs::storage::exv::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(i))).evaluate(0);
            double fjx = form_factor::foxs::storage::exv::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(j))).evaluate(0);

            count.at( std::round(aa_distances(i, j)*2)) += 2;
            xx_sum.at(std::round(aa_distances(i, j)*2)) += 2*fix*fjx;
            aa_sum.at(std::round(aa_distances(i, j)*2)) += 2*fi*fj;
        }
    }
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        count.at(0) += 1;
        xx_sum.at(0) += std::pow(form_factor::foxs::storage::exv::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(i))).evaluate(0), 2);
        aa_sum.at(0) += std::pow(form_factor::foxs::storage::atomic::get_form_factor(static_cast<form_factor::form_factor_t>(cc.get_ff_type(i))).evaluate(0), 2);
    }

    for (unsigned int i = 0; i < 10; ++i) {
        std::cout << "d = " << constants::axes::d_axis.get_bin_value(i) << std::endl;
        std::cout << "\taa_sum[" << i << "] = " << aa_sum[i] << std::endl;
        std::cout << "\txx_sum[" << i << "] = " << xx_sum[i] << std::endl;
        std::cout << "\t count[" << i << "] = " << count[i] << std::endl;
    }
    std::cout << "Exiting at DebugManager.cpp:debye_transform()" << std::endl;
    exit(1);

    for (const auto& q : q_axis) {
        double Iq = 0;
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                const auto& atom_i = atoms[i];
                const auto& atom_j = atoms[j];

                // atom
                double Zi = atom_i.get_effective_charge()*atom_i.get_occupancy();
                double Zj = atom_j.get_effective_charge()*atom_j.get_occupancy();
                double fi = Zi*form_factor::storage::atomic::get_form_factor(form_factor::get_type(atom_i.get_element(), atom_i.get_atomic_group())).evaluate(q);
                double fj = Zj*form_factor::storage::atomic::get_form_factor(form_factor::get_type(atom_j.get_element(), atom_i.get_atomic_group())).evaluate(q);
                double qr = q*aa_distances(i, j);

                // exv
                double fix = form_factor::storage::exv::get_form_factor(form_factor::get_type(atom_i.get_element(), atom_i.get_atomic_group())).evaluate(q);
                double fjx = form_factor::storage::exv::get_form_factor(form_factor::get_type(atom_j.get_element(), atom_i.get_atomic_group())).evaluate(q);

                double tmp = Zi*Zj*fi*fj + fix*fjx - 2*Zi*fi*fix;
                if (qr < 1e-9) {
                    Iq += tmp;
                } else {
                    Iq += tmp*sin(qr)/qr;
                }
            }

            for (unsigned int j = 0; j < protein->water_size(); ++j) {
                const auto& atom_i = atoms[i];
                const auto& water = protein->get_waters()[j];

                // water
                double Zi = atom_i.get_effective_charge()*atom_i.get_occupancy();
                double Zj = water.get_effective_charge()*water.get_occupancy()*cw;
                double fi = Zi*form_factor::storage::atomic::get_form_factor(form_factor::get_type(atom_i.get_element(), atom_i.get_atomic_group())).evaluate(q);
                double fj = Zj*form_factor::storage::atomic::O.evaluate(q);
                double qr = q*aw_distances(i, j);

                // exv
                double Z_exv_i = Z_exv*atom_i.get_occupancy()*cx;
                double f_exv = Z_exv_i*form_factor::storage::atomic::excluded_volume.evaluate(q);

                double tmp = 2*Zi*Zj*fi*fj - 2*Zi*Z_exv_i*fi*f_exv;
                if (qr < 1e-9) {
                    Iq += tmp;
                } else {
                    Iq += tmp*sin(qr)/qr;
                }
            }
        }

        for (unsigned int i = 0; i < protein->water_size(); ++i) {
            for (unsigned int j = 0; j < protein->water_size(); ++j) {
                const auto& water_i = protein->get_waters()[i];
                const auto& water_j = protein->get_waters()[j];

                double Zi = water_i.get_effective_charge()*water_i.get_occupancy()*cw;
                double Zj = water_j.get_effective_charge()*water_j.get_occupancy()*cw;
                double fw = Zj*form_factor::storage::atomic::O.evaluate(q);
                double qr = q*ww_distances(i, j);

                double tmp = Zi*Zj*fw*fw;
                if (qr < 1e-9) {
                    Iq += tmp;
                } else {
                    Iq += tmp*sin(qr)/qr;
                }
            }
        }

        I.push_back(Iq);
    }
    return ScatteringProfile(I, Axis(settings::axes::qmin, settings::axes::qmax, 100));
}

DebugManager::~DebugManager() = default;

std::unique_ptr<DistanceHistogram> DebugManager::calculate() {
    return std::make_unique<DebugDistanceHistogram>(protein);
}

std::unique_ptr<CompositeDistanceHistogram> DebugManager::calculate_all() {
    return std::make_unique<DebugDistanceHistogram>(protein);
}