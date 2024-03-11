#include "hist/Histogram.h"
#include "plots/PlotDataset.h"
#include "plots/PlotHistogram.h"
#include "settings/EMSettings.h"
#include <CLI/CLI.hpp>

#include <em/ImageStack.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/File.h>
#include <rigidbody/sequencer/Sequencer.h>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <data/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>
#include <settings/MoleculeSettings.h>
#include <fitter/HydrationFitter.h>
#include <form_factor/FormFactor.h>
#include <hydrate/Grid.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <em/detail/header/MRCHeader.h>
#include <em/detail/header/data/MRCData.h>

// fix test
int main(int argc, char const *argv[]) {
    settings::general::verbose = true;

    data::Molecule protein_2epe("test/files/2epe.pdb");
    data::Molecule protein_LAR12("test/files/LAR1-2.pdb");
    protein_2epe.generate_new_hydration();
    protein_LAR12.generate_new_hydration();

    fitter::HydrationFitter fitter("test/files/2epe.dat", protein_2epe.get_histogram());
    double chi2 = fitter.fit()->fval;

    std::cout << "chi2: " << chi2 << std::endl;

    fitter.set_scattering_hist(protein_LAR12.get_total_histogram());
    double _chi2 = fitter.fit()->fval;

    std::cout << "chi2 " << _chi2 << std::endl;

    fitter.set_scattering_hist(protein_2epe.get_total_histogram());
    _chi2 = fitter.fit()->fval;

    std::cout << "chi2 " << _chi2 << std::endl;
}

// generate histogram over charge values in EM map
// int main(int argc, char const *argv[]) {
//     std::string map = argv[1];

//     em::ImageStack images(map);
//     std::vector<double> bins(1e5);
//     Axis x(-2, 10, 1e5);
//     auto inv_rms = 1./images.rms();
//     for (const auto& image : images.images()) {
//         for (auto val : image.get_data()) {
//             int index = std::floor((val*inv_rms - x.min)/x.width());
//             if (index < 0 || bins.size() < index) {
//                 continue;
//             }
//             bins[index]++;
//         }
//     }

//     hist::Histogram dataset(std::move(bins), x);
//     plots::PlotHistogram::quick_plot(dataset, plots::PlotOptions(style::draw::points, {{"logy", true}}), "output/stuff/charge_hist.png");
// }

// int main(int argc, char const *argv[]) {
//     settings::molecule::center = false;
//     settings::molecule::use_effective_charge = false; 
//     settings::general::threads = 6;
//     settings::em::hydrate = false;
//     settings::em::fixed_weights = true;

//     // generate big sphere
//     auto lims = Limit3D(-50, 50, -50, 50, -50, 50);
//     grid::Grid grid(lims);
//     double radius = 25;
//     double radius2 = radius*radius;
//     auto axes = grid.get_axes();
//     Vector3<double> center = grid.to_xyz(grid.get_center());
//     std::cout << center << std::endl;
//     for (unsigned int i = 0; i < axes.x.bins; ++i) {
//         for (unsigned int j = 0; j < axes.y.bins; ++j) {
//             for (unsigned int k = 0; k < axes.z.bins; ++k) {
//                 if (grid.to_xyz(i, j, k).distance2(center) < radius2) {
//                     grid.grid.index(i, j, k) = grid::detail::VOLUME;
//                 }
//             }
//         }
//     }
//     auto loc = "temp/test/em/sphere.pdb";
//     grid.save(loc);

//     data::Molecule protein(loc);
//     auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
//     Iq.as_dataset().save("temp/test/em/sphere_Iq.dat");

//     std::unique_ptr<em::detail::header::MRCHeader> header;
//     {
//         em::detail::header::MRCData header_data;
//         header_data.cella_x = axes.x.span();
//         header_data.cella_y = axes.y.span();
//         header_data.cella_z = axes.z.span();
//         header_data.nx = axes.x.bins;
//         header_data.ny = axes.y.bins;
//         header_data.nz = axes.z.bins;
//         header = std::make_unique<em::detail::header::MRCHeader>(std::move(header_data));
//     }

//     std::vector<em::Image> images(lims.z.span()/settings::grid::width, Matrix<float>(0, 0));
//     for (unsigned int k = 0; k < images.size(); ++k) {
//         Matrix<float> data(axes.x.bins, axes.y.bins);
//         for (unsigned int i = 0; i < axes.x.bins; ++i) {
//             for (unsigned int j = 0; j < axes.y.bins; ++j) {
//                 double dist = std::sqrt(grid.to_xyz(i, j, k).distance2(center));
//                 data.index(i, j) = radius/dist;
//             }
//         }
//         data.index(50, 50) = 2*data.index(50, 51);
//         images[k] = em::Image(data, header.get(), k);
//     }

//     em::ImageStack stack(images);
//     // auto[fit, landscape] = stack.cutoff_scan_fit(100, std::move(Iq));
//     auto fit = stack.fit("temp/test/em/sphere_Iq.dat");
//     // plots::PlotLandscape::quick_plot(landscape, "temp/test/em/sphere_landscape.png");
//     std::cout << fit->fval/fit->dof << std::endl;
// }

// int main(int argc, char const *argv[]) {
//     auto exact = [] (const data::Molecule& molecule, double exv_radius) {
//         container::Container2D<double> distances(molecule.get_atoms().size(), molecule.get_atoms().size());
//         auto atoms = molecule.get_atoms();
//         for (unsigned int i = 0; i < atoms.size(); ++i) {
//             for (unsigned int j = 0; j < atoms.size(); ++j) {
//                 distances(i, j) = atoms[i].distance(atoms[j]);
//             }
//         }

//         auto qaxis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
//         auto q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
//         form_factor::FormFactor ff = form_factor::ExvFormFactor(std::pow(2*exv_radius, 3));
//         hist::ScatteringProfile I(qaxis);
//         for (unsigned int q = q0; q < q0+qaxis.bins; ++q) {
//             double sum = 0;
//             for (unsigned int i = 0; i < atoms.size(); ++i) {
//                 for (unsigned int j = 0; j < atoms.size(); ++j) {
//                     double qd = constants::axes::q_vals[q]*distances(i, j);
//                     if (qd < 1e-6) {
//                         sum += std::pow(ff.evaluate(constants::axes::q_vals[q]), 2);
//                     } else {
//                         sum += std::pow(ff.evaluate(constants::axes::q_vals[q]), 2)*std::sin(qd)/qd;
//                     }
//                 }
//             }
//             I.index(q-q0) = sum;
//         }
//         return I;
//     };

//     settings::axes::qmin = 5e-2; 
//     settings::axes::qmax = 1;
//     settings::grid::width = 1;
//     settings::grid::exv_radius = 1;
//     settings::hist::weighted_bins = true;
//     settings::general::output = "temp/stuff/comparison/";
//     data::Molecule protein("data/6lyz/6lyz.pdb");
//     std::vector<SimpleDataset> profiles;
//     for (double rx = 0.5; rx <= 3; rx += 0.5) {
//         settings::grid::exv_radius = rx;
//         protein.clear_grid();
//         hist::CompositeDistanceHistogramFFGrid::regenerate_table();
//         auto h = hist::HistogramManagerMTFFGrid<true>(&protein).calculate_all();
//         auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
//         auto profile = h_cast->get_profile_xx();
//         profiles.push_back(profile.as_dataset());
//     }

//     auto jan_1 = SimpleDataset("temp/stuff/comparison/jan_1.dat");
//     auto jan_2 = SimpleDataset("temp/stuff/comparison/jan_2.dat");
//     auto jan_3 = SimpleDataset("temp/stuff/comparison/jan_3.dat");

//     profiles[0].normalize(1);
//     profiles[1].normalize(10);
//     profiles[2].normalize(100);
//     jan_1.normalize(1);
//     jan_2.normalize(10);
//     jan_3.normalize(100);

//     plots::PlotIntensity()
//         .plot(profiles[0], plots::PlotOptions({{"legend", "1"},     {"lw", 2}, {"color", style::color::red}, {"ylimits", Limit{1e-4, 110}}}))
//         .plot(jan_1,       plots::PlotOptions({{"legend", "Jan 1"}, {"lw", 2}, {"color", style::color::red}, {"linestyle", style::line::dashed}}))
//         .plot(profiles[1], plots::PlotOptions({{"legend", "2"},     {"lw", 2}, {"color", style::color::blue}}))
//         .plot(jan_2,       plots::PlotOptions({{"legend", "Jan 2"}, {"lw", 2}, {"color", style::color::blue}, {"linestyle", style::line::dashed}}))
//         .plot(profiles[2], plots::PlotOptions({{"legend", "3"},     {"lw", 2}, {"color", style::color::green}}))
//         .plot(jan_3,       plots::PlotOptions({{"legend", "Jan 3"}, {"lw", 2}, {"color", style::color::green}, {"linestyle", style::line::dashed}}))
//     .save("temp/stuff/comparison/compare.png");

//     constants::radius::set_dummy_radius(0);
//     settings::grid::rvol = 0;

//     settings::grid::width = 1;
//     settings::grid::exv_radius = 0.5;
//     settings::grid::save_exv = true;
//     hist::CompositeDistanceHistogramFFGrid::regenerate_table();
//     data::Molecule jan1("temp/stuff/comparison/jan_1.pdb");
//     for (auto& b : jan1.get_bodies()) {for (auto& a : b.get_atoms()){a.element = constants::atom_t::dummy;}}
//     auto p1 = static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<true>(&jan1).calculate_all().get())->get_profile_xx().as_dataset();
//     // auto exact1 = exact(jan1, 0.5).as_dataset();
//     settings::grid::save_exv = false;

//     constants::radius::set_dummy_radius(1);
//     settings::grid::rvol = 1;

//     settings::grid::exv_radius = 1;
//     hist::CompositeDistanceHistogramFFGrid::regenerate_table();
//     data::Molecule jan2("temp/stuff/comparison/jan_2.pdb");
//     for (auto& b : jan2.get_bodies()) {for (auto& a : b.get_atoms()){a.element = constants::atom_t::dummy;}}
//     auto p2 = static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<true>(&jan2).calculate_all().get())->get_profile_xx().as_dataset();
//     auto exact2 = exact(jan2, 1).as_dataset();

//     constants::radius::set_dummy_radius(2);
//     settings::grid::rvol = 2;

//     settings::grid::exv_radius = 1.5;
//     hist::CompositeDistanceHistogramFFGrid::regenerate_table();
//     data::Molecule jan3("temp/stuff/comparison/jan_3.pdb");
//     for (auto& b : jan3.get_bodies()) {for (auto& a : b.get_atoms()){a.element = constants::atom_t::dummy;}}
//     auto p3 = static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<true>(&jan3).calculate_all().get())->get_profile_xx().as_dataset();
//     auto exact3 = exact(jan3, 1.5).as_dataset();

//     p1.normalize(1);
//     // exact1.normalize(1);
//     p2.normalize(10);
//     exact2.normalize(10);
//     p3.normalize(100);
//     exact3.normalize(100);

//     plots::PlotIntensity()
//         .plot(p1,     plots::PlotOptions({{"legend", "Grid 1"}, {"lw", 2}, {"color", style::color::red}, {"ylimits", Limit{1e-4, 110}}}))
//         .plot(jan_1,  plots::PlotOptions({{"legend", "Jan 1"}, {"lw", 2}, {"color", style::color::red}, {"linestyle", style::line::dashed}}))
//         // .plot(exact1, plots::PlotOptions({{"legend", "Exact 1"}, {"lw", 2}, {"color", style::color::red}, {"linestyle", style::line::dotted}}))
//         .plot(p2,    plots::PlotOptions({{"legend", "Grid 2"}, {"lw", 2}, {"color", style::color::blue}}))
//         .plot(jan_2, plots::PlotOptions({{"legend", "Jan 2"}, {"lw", 2}, {"color", style::color::blue}, {"linestyle", style::line::dashed}}))
//         .plot(exact2, plots::PlotOptions({{"legend", "Exact 2"}, {"lw", 2}, {"color", style::color::blue}, {"linestyle", style::line::dotted}}))
//         .plot(p3,    plots::PlotOptions({{"legend", "Grid 3"}, {"lw", 2}, {"color", style::color::green}}))
//         .plot(jan_3, plots::PlotOptions({{"legend", "Jan 3"}, {"lw", 2}, {"color", style::color::green}, {"linestyle", style::line::dashed}}))
//         .plot(exact3, plots::PlotOptions({{"legend", "Exact 3"}, {"lw", 2}, {"color", style::color::green}, {"linestyle", style::line::dotted}}))
//     .save("temp/stuff/comparison/jan.png");
// }

// int main(int argc, char const *argv[]) {
//     settings::axes::qmin = 5e-2; settings::axes::qmax = 1;
//     settings::grid::width = 1;
//     settings::grid::exv_radius = 1;
//     data::Molecule protein("data/6lyz/6lyz.pdb");
//     plots::PlotIntensity plot;
//     double c = 0;
//     for (double rx = 0.5; rx <= 1.5; rx += 0.5) {
//         settings::grid::exv_radius = rx;
//         protein.clear_grid();
//         hist::CompositeDistanceHistogramFFGrid::regenerate_table();
//         auto h = hist::HistogramManagerMTFFGrid<true>(&protein).calculate_all();
//         auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
//         auto profile = h_cast->get_profile_xx().as_dataset();
//         if (c == 0) {
//             c = profile.normalize();
//         } else {
//             profile.scale_y(c);
//         }
//         plot.plot(profile, plots::PlotOptions({{"legend", utility::round(rx, 2)}, {"lw", 2}, {"color", style::color::next()}, {"yrange", Limit(1e-4, 1.1)}}));
//     }
//     plot.save("temp/stuff/xx_variation.png");
// }

//************************************************************************************************* 
//********************************* PLOT SCATTERING STUFF *****************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::axes::qmin = 0.02;
//     settings::axes::qmax = 1;

//     std::string base_path = "temp/debug/";
//     SimpleDataset I_aa(base_path + "ausaxs/ausaxs_aa.dat");
//     SimpleDataset I_ax(base_path + "ausaxs/ausaxs_ax.dat");
//     SimpleDataset I_xx(base_path + "ausaxs/ausaxs_xx.dat");
//     SimpleDataset I_aw(base_path + "ausaxs/ausaxs_aw.dat");
//     SimpleDataset I_wx(base_path + "ausaxs/ausaxs_wx.dat");
//     SimpleDataset I_ww(base_path + "ausaxs/ausaxs_ww.dat");

//     SimpleDataset coords(base_path + "COORDINATE.dat");
//     SimpleDataset C_xx(base_path + "exclvol.dat");
//     SimpleDataset C_aa(base_path + "vaccum.dat");
//     SimpleDataset foxs(base_path + "foxs_vacuum.dat");

//     SimpleDataset foxs_aa(base_path + "foxs_aa.dat");
//     SimpleDataset foxs_ax(base_path + "foxs_ax.dat");
//     SimpleDataset foxs_xx(base_path + "foxs_xx.dat");
//     SimpleDataset foxs_aw(base_path + "foxs_aw.dat");
//     SimpleDataset foxs_wx(base_path + "foxs_wx.dat");
//     SimpleDataset foxs_ww(base_path + "foxs_ww.dat");
    
//     // SimpleDataset lyshr(base_path + "LYSHR.RSR");

//     I_aa.normalize(1);
//     I_ax.normalize(1);
//     I_xx.normalize(1);
//     coords.normalize(1);
//     C_xx.normalize(1);
//     C_aa.normalize(1);
//     foxs.normalize(1);

//     SimpleDataset I_sum(I_aa);
//     for (unsigned int i = 0; i < I_aa.size(); ++i) {
//         I_aa.y(i) = std::abs(I_aa.y(i));
//         I_ax.y(i) = std::abs(I_ax.y(i));
//         I_xx.y(i) = std::abs(I_xx.y(i));
//         I_aw.y(i) = std::abs(I_aw.y(i));
//         I_wx.y(i) = std::abs(I_wx.y(i));
//         I_ww.y(i) = std::abs(I_ww.y(i));

//         double cy = coords.interpolate_y(I_aa.x(i));
//         I_sum.y(i) = I_aa.y(i) + I_ax.y(i) + I_xx.y(i);
//         // I_aa.y(i) /= cy;
//         // I_ax.y(i) /= cy;
//         // I_xx.y(i) /= cy;
//     }

//     for (unsigned int i = 0; i < foxs_aa.size(); ++i) {
//         foxs_aa.y(i) = std::abs(foxs_aa.y(i));
//         foxs_ax.y(i) = std::abs(foxs_ax.y(i));
//         foxs_xx.y(i) = std::abs(foxs_xx.y(i));
//         foxs_aw.y(i) = std::abs(foxs_aw.y(i));
//         foxs_wx.y(i) = std::abs(foxs_wx.y(i));
//         foxs_ww.y(i) = std::abs(foxs_ww.y(i));
//     }

//     I_aa.normalize(1);
//     I_ax.normalize(1);
//     I_xx.normalize(1);
//     I_aw.normalize(1);
//     I_wx.normalize(1);
//     I_ww.normalize(1);

//     foxs_aa.normalize(1);
//     foxs_ax.normalize(1);
//     foxs_xx.normalize(1);
//     foxs_aw.normalize(1);
//     foxs_wx.normalize(1);
//     foxs_ww.normalize(1);

//     coords.normalize(1);
//     C_xx.normalize(1);
//     C_aa.normalize(1);
//     I_sum.normalize(1);
//     foxs.normalize(1);

//     I_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}});
//     I_ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ax"}, {"color", style::color::blue}});
//     I_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::green}});
//     I_aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aw"}, {"color", style::color::pink}});
//     I_wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xw"}, {"color", style::color::purple}});
//     I_ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ww"}, {"color", style::color::brown}});

//     foxs_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aa"}, {"color", style::color::red}});
//     foxs_ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ax"}, {"color", style::color::blue}});
//     foxs_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xx"}, {"color", style::color::green}});
//     foxs_aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aw"}, {"color", style::color::pink}});
//     foxs_wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xw"}, {"color", style::color::purple}});
//     foxs_ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ww"}, {"color", style::color::brown}});

//     coords.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "coords"}, {"linestyle", style::line::dashed}, {"color", style::color::red}});
//     C_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_xx"}, {"linestyle", style::line::dashed}, {"color", style::color::blue}});
//     C_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_aa"}, {"linestyle", style::line::dashed}, {"color", style::color::green}});
//     I_sum.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "I_sum"}, {"color", style::color::black}});
//     foxs.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "foxs"}, {"color", style::color::black}});
//     // lyshr.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "lyshr"}, {"color", style::color::black}});

//     // I_aa.save(base_path + "AUSAXS_aa.dat");
//     // I_ax.save(base_path + "AUSAXS_ax.dat");
//     // I_xx.save(base_path + "AUSAXS_xx.dat");
//     // exclvol.save(base_path + "scaled_exclvol.dat");
//     // vacuum.save(base_path + "scaled_vacuum.dat");

//     plots::PlotDataset(I_aa)
//         .plot(I_xx)
//         .plot(C_xx)
//         .plot(C_aa)
//     .save(base_path + "compare.png");

//     plots::PlotDataset(I_aa)
//         .plot(C_aa)
//         .plot(foxs)
//     .save(base_path + "vacuum.png");

//     plots::PlotDataset(I_aa)
//         .plot(I_xx)
//         .plot(I_ww)
//         .plot(foxs_aa)
//         .plot(foxs_xx)
//         .plot(foxs_ww)
//     .save(base_path + "compare_foxs.png");

//     plots::PlotDataset(I_ax)
//         .plot(I_aw)
//         .plot(I_wx)
//         .plot(foxs_ax)
//         .plot(foxs_aw)
//         .plot(foxs_wx)
//     .save(base_path + "compare_foxs_cross.png");

//     plots::PlotDataset(foxs_xx)
//         .plot(C_xx)
//         .plot(I_xx)
//     .save(base_path + "exv.png");

//     for (unsigned int i = 0; i < C_xx.size(); ++i) {
//         C_aa.y(i) /= I_aa.interpolate_y(C_aa.x(i));
//         C_xx.y(i) /= I_xx.interpolate_y(C_xx.x(i));
//     }
//     // C_xx.normalize();
//     // C_aa.normalize();

//     for (unsigned int i = 0; i < foxs_xx.size(); ++i) {
//         foxs_aa.y(i) /= I_aa.interpolate_y(foxs_aa.x(i));
//         foxs_xx.y(i) /= I_xx.interpolate_y(foxs_xx.x(i));
//     }
//     // foxs_aa.normalize();
//     // foxs_xx.normalize();

//     foxs_xx.add_plot_options({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}});
//     plots::PlotDataset(foxs_xx)
//         .plot(C_xx)
//         .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
//     .save(base_path + "exv_normalized.png");

//     foxs_aa.add_plot_options({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}});
//     plots::PlotDataset(foxs_aa)
//         .plot(C_aa)
//         .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
//     .save(base_path + "vacuum_normalized.png");


//     // for (unsigned int i = 0; i < foxs_xx.size(); ++i) {foxs_xx.y(i) /= foxs_aa.y(i);}
//     // for (unsigned int i = 0; i < I_xx.size(); ++i) {I_xx.y(i) /= I_aa.y(i);}

//     // plots::PlotDataset(foxs_xx)
//     //     .plot(I_xx)
//     //     .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
//     // .save(base_path + "xx_div_aa.png");
// }

//*************************************************************************************************
//********************************* PLOT FORM FACTORS *********************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::general::output = "temp/stuff/ff/";
//     auto plot = [] (form_factor::FormFactor ff, std::string name) {
//         auto q = Axis(0, 1, 100).as_vector();
//         std::vector<double> y(q.size());
//         for (unsigned int i = 0; i < q.size(); ++i) {
//             y[i] = ff.evaluate(q[i]);
//         }
//         SimpleDataset d(q, y);
//         d.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff"}});
//         return d;
//     };

//     auto C = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::C), "carbon");
//     auto O = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::O), "oxygen");
//     auto EV = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::EXCLUDED_VOLUME), "excluded_volume");
//     auto other = plot(form_factor::storage::get_form_factor(form_factor::form_factor_t::OTHER), "other");
//     C.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff"}, {"legend", "Carbon"}, {"color", "red"}});
//     O.add_plot_options({{"legend", "Oxygen"}, {"color", "blue"}});
//     EV.add_plot_options({{"legend", "Excluded volume"}, {"color", "green"}});
//     other.add_plot_options({{"legend", "Other"}, {"color", "black"}});

//     plots::PlotDataset(C)
//         .plot(O)
//         .plot(EV)
//         .plot(other)
//     .save(settings::general::output + "ff.png");
// }

// int main(int argc, char const *argv[]) {
//     settings::protein::use_effective_charge = false;
//     io::ExistingFile mapfile(argv[1]);
//     io::ExistingFile mfile(argv[2]);
//     em::ImageStack images(mapfile);

//     settings::em::alpha_levels = {14, 15};
//     settings::fit::max_iterations = 100;
//     Axis axis(settings::em::alpha_levels, settings::fit::max_iterations);
//     fitter::HydrationFitter fitter(mfile);
//     unsigned int c = 0;
//     for (auto& level : axis.as_vector()) {
//         auto protein = images.get_protein(level);
//         protein->generate_new_hydration();
//         fitter.set_scattering_hist(protein->get_histogram());
//         auto res = fitter.fit();
//         std::cout << "Step " << ++c << ": Evaluated cutoff value " << level << " with chi2 " << res->fval << std::endl;
//     }
// }

//*************************************************************************************************
//********************************* PLOT EXCLUDED VOLUME ******************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::protein::use_effective_charge = false;
//     io::ExistingFile file(argv[1]);
//     Protein protein(file);
//     protein.get_grid()->expand_volume();
//     protein.get_grid()->save("temp/stuff/grid.pdb");
// }

// int main(int argc, char const *argv[]) {
//     settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
//     settings::protein::use_effective_charge = false;

//         int distance = settings::rigidbody::bond_distance;
//         Atom a1 = Atom(Vector3<double>(0, 0, 0*distance), 1, "C", "C", 1);
//         Atom a2 = Atom(Vector3<double>(0, 0, 1*distance), 1, "C", "C", 1);
//         Atom a3 = Atom(Vector3<double>(0, 0, 2*distance), 1, "C", "C", 1);
//         Atom a4 = Atom(Vector3<double>(0, 0, 3*distance), 1, "C", "C", 1);

//         Body b1 = Body(std::vector<Atom>{a1});
//         Body b2 = Body(std::vector<Atom>{a2});
//         Body b3 = Body(std::vector<Atom>{a3});
//         Body b4 = Body(std::vector<Atom>{a4});
//         std::vector<Body> ap = {b1, b2, b3, b4};
//         rigidbody::RigidBody rigidbody(ap);
//         assert(rigidbody.get_constraint_manager()->distance_constraints.size() == 3);

//         rigidbody::RigidBody rigidbody2 = rigidbody::BodySplitter::split("data/LAR1-2/LAR1-2.pdb", {9, 99});
//         assert(rigidbody2.get_constraint_manager()->distance_constraints.size() == 2);
// }

//*************************************************************************************************
//********************************* DEBUG SEQUENCER ***********************************************
//*************************************************************************************************
// int main(int argc, char const *argv[]) {
//     settings::grid::scaling = 2;
//     settings::grid::cubic = true;
//     settings::general::verbose = true;
//     settings::axes::distance_bin_width = 0.1;

//     CLI::App app{"Rigid-body optimization."};
//     io::File pdb, mfile, settings;
//     std::vector<unsigned int> constraints;
//     app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
//     app.add_option("input_m", mfile, "Path to the measuremed data.")->required()->check(CLI::ExistingFile);
//     app.add_option("output", settings::general::output, "Path to save the hydrated file at.")->default_val("output/rigidbody/");
//     auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile);
//     app.add_option("--reduce,-r", settings::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
//     app.add_option("--grid_width,-w", settings::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
//     app.add_option("--bin_width", settings::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
//     app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
//     app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
//     auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
//     app.add_option("--iterations", settings::rigidbody::iterations, "Maximum number of iterations. Default: 1000.");
//     app.add_option("--constraints", settings::rigidbody::detail::constraints, "Constraints to apply to the rigid body.");
//     app.add_flag("--center,!--no-center", settings::protein::center, "Decides whether the protein will be centered. Default: true.");
//     app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
//     CLI11_PARSE(app, argc, argv);
    
//     //###################//
//     //### PARSE INPUT ###//
//     //###################//
//     settings::general::output += mfile.stem() + "/";

//     // if a settings file was provided
//     if (p_settings->count() != 0) {
//         settings::read(settings);        // read it
//         CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
//     } else {                            // otherwise check if there is a settings file in the same directory
//         if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
//             CLI11_PARSE(app, argc, argv);
//         }
//     }
//     if (settings::rigidbody::detail::constraints.empty()) {
//         throw except::missing_option("rigidbody: Constraints must be supplied. Use --constraints to specify them.");
//     }

//     rigidbody::sequencer::Sequencer(mfile, rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints))
//         .body_select_strategy(settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect)
//         .parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly)
//             .decay_strategy(settings::rigidbody::DecayStrategyChoice::Exponential)
//             .amplitude(0.5)
//         .transform_strategy(settings::rigidbody::TransformationStrategyChoice::RigidTransform)
//         .loop(5)
//     .execute();

//     return 0;
// }