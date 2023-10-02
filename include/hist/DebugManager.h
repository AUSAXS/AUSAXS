#pragma once

#include <hist/HistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/detail/FormFactor.h>
#include <hist/DebyeLookupTable.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/Protein.h>
#include <settings/All.h>
#include <utility/Container2D.h>

namespace hist {
    class DebugDistanceHistogram : public CompositeDistanceHistogramFF {
        public: 
            DebugDistanceHistogram(Protein* protein) : protein(protein) {
                settings::axes::bins = 100;
                this->q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
            }

            ScatteringHistogram debye_transform() const override {
                static unsigned int counter = 0;
                std::cout << "debye transform " << counter++ << std::endl;

                std::vector<Atom> atoms = protein->get_atoms();
                double Z_exv = protein->get_excluded_volume()*constants::charge::density::water/atoms.size();

                std::vector<double> I;
                I.reserve(q_axis.size());

                Container2D<double> aa_distances(atoms.size(), atoms.size());
                Container2D<double> aw_distances(atoms.size(), protein->get_waters().size());
                Container2D<double> ww_distances(protein->get_waters().size(), protein->get_waters().size());
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

                for (const auto& q : q_axis) {
                    double Iq = 0;
                    for (unsigned int i = 0; i < atoms.size(); ++i) {
                        for (unsigned int j = 0; j < atoms.size(); ++j) {
                            const auto& atom_i = atoms[i];
                            const auto& atom_j = atoms[j];

                            // atom
                            double Zi = atom_i.get_effective_charge()*atom_i.get_occupancy();
                            double Zj = atom_j.get_effective_charge()*atom_j.get_occupancy();
                            double fi = Zi*hist::detail::FormFactorStorage::get_form_factor(hist::detail::FormFactor::get_type(atom_i)).evaluate(q);
                            double fj = Zj*hist::detail::FormFactorStorage::get_form_factor(hist::detail::FormFactor::get_type(atom_j)).evaluate(q);
                            double qr = q*aa_distances(i, j);

                            // exv
                            double Z_exv_i = Z_exv*atom_i.get_occupancy()*exv_scaling;
                            double Z_exv_j = Z_exv*atom_j.get_occupancy()*exv_scaling;
                            double f_exv = Z_exv_i*hist::detail::FormFactorStorage::excluded_volume.evaluate(q);

                            double tmp = Zi*Zj*fi*fj + Z_exv_i*Z_exv_j*f_exv*f_exv - 2*Zi*Z_exv_i*fi*f_exv;
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
                            double Zj = water.get_effective_charge()*water.get_occupancy()*w_scaling;
                            double fi = Zi*hist::detail::FormFactorStorage::get_form_factor(hist::detail::FormFactor::get_type(atom_i)).evaluate(q);
                            double fj = Zj*hist::detail::FormFactorStorage::oxygen.evaluate(q);
                            double qr = q*aw_distances(i, j);

                            // exv
                            double Z_exv_i = Z_exv*atom_i.get_occupancy()*exv_scaling;
                            double f_exv = Z_exv_i*hist::detail::FormFactorStorage::excluded_volume.evaluate(q);

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

                            double Zi = water_i.get_effective_charge()*water_i.get_occupancy()*w_scaling;
                            double Zj = water_j.get_effective_charge()*water_j.get_occupancy()*w_scaling;
                            double fw = Zj*hist::detail::FormFactorStorage::oxygen.evaluate(q);
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
                return ScatteringHistogram(I, Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins));
            }

        private:
            Protein* protein;
    };

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class DebugManager : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			~DebugManager() override = default;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override {
                return std::make_unique<DebugDistanceHistogram>(protein);
            }

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<CompositeDistanceHistogram> calculate_all() override {
                return std::make_unique<DebugDistanceHistogram>(protein);
            }
	};
}