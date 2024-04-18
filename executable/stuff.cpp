#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

int main(int argc, char const *argv[]) {
    settings::axes::qmin = 0.02;
    settings::axes::qmax = 1;

    std::string base_path = "temp/debug/";
    SimpleDataset I_aa(base_path + "ausaxs/ausaxs_aa.dat");
    SimpleDataset I_ax(base_path + "ausaxs/ausaxs_ax.dat");
    SimpleDataset I_xx(base_path + "ausaxs/ausaxs_xx.dat");
    SimpleDataset I_aw(base_path + "ausaxs/ausaxs_aw.dat");
    SimpleDataset I_wx(base_path + "ausaxs/ausaxs_wx.dat");
    SimpleDataset I_ww(base_path + "ausaxs/ausaxs_ww.dat");

    SimpleDataset coords(base_path + "COORDINATE.dat");
    SimpleDataset C_xx(base_path + "exclvol.dat");
    SimpleDataset C_aa(base_path + "vaccum.dat");
    SimpleDataset foxs(base_path + "foxs_vacuum.dat");

    SimpleDataset foxs_aa(base_path + "foxs_aa.dat");
    SimpleDataset foxs_ax(base_path + "foxs_ax.dat");
    SimpleDataset foxs_xx(base_path + "foxs_xx.dat");
    SimpleDataset foxs_aw(base_path + "foxs_aw.dat");
    SimpleDataset foxs_wx(base_path + "foxs_wx.dat");
    SimpleDataset foxs_ww(base_path + "foxs_ww.dat");
    
    // SimpleDataset lyshr(base_path + "LYSHR.RSR");

    I_aa.normalize(1);
    I_ax.normalize(1);
    I_xx.normalize(1);
    coords.normalize(1);
    C_xx.normalize(1);
    C_aa.normalize(1);
    foxs.normalize(1);

    SimpleDataset I_sum(I_aa);
    for (unsigned int i = 0; i < I_aa.size(); ++i) {
        I_aa.y(i) = std::abs(I_aa.y(i));
        I_ax.y(i) = std::abs(I_ax.y(i));
        I_xx.y(i) = std::abs(I_xx.y(i));
        I_aw.y(i) = std::abs(I_aw.y(i));
        I_wx.y(i) = std::abs(I_wx.y(i));
        I_ww.y(i) = std::abs(I_ww.y(i));

        double cy = coords.interpolate_x(I_aa.x(i), 1);
        I_sum.y(i) = I_aa.y(i) + I_ax.y(i) + I_xx.y(i);
        // I_aa.y(i) /= cy;
        // I_ax.y(i) /= cy;
        // I_xx.y(i) /= cy;
    }

    for (unsigned int i = 0; i < foxs_aa.size(); ++i) {
        foxs_aa.y(i) = std::abs(foxs_aa.y(i));
        foxs_ax.y(i) = std::abs(foxs_ax.y(i));
        foxs_xx.y(i) = std::abs(foxs_xx.y(i));
        foxs_aw.y(i) = std::abs(foxs_aw.y(i));
        foxs_wx.y(i) = std::abs(foxs_wx.y(i));
        foxs_ww.y(i) = std::abs(foxs_ww.y(i));
    }

    I_aa.normalize(1);
    I_ax.normalize(1);
    I_xx.normalize(1);
    I_aw.normalize(1);
    I_wx.normalize(1);
    I_ww.normalize(1);

    foxs_aa.normalize(1);
    foxs_ax.normalize(1);
    foxs_xx.normalize(1);
    foxs_aw.normalize(1);
    foxs_wx.normalize(1);
    foxs_ww.normalize(1);

    coords.normalize(1);
    C_xx.normalize(1);
    C_aa.normalize(1);
    I_sum.normalize(1);
    foxs.normalize(1);

    // I_aa.save(base_path + "AUSAXS_aa.dat");
    // I_ax.save(base_path + "AUSAXS_ax.dat");
    // I_xx.save(base_path + "AUSAXS_xx.dat");
    // exclvol.save(base_path + "scaled_exclvol.dat");
    // vacuum.save(base_path + "scaled_vacuum.dat");

    plots::PlotDataset(I_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}}))
        .plot(I_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::green}}))
        .plot(C_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_xx"}, {"linestyle", style::line::dashed}, {"color", style::color::blue}}))
        .plot(C_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_aa"}, {"linestyle", style::line::dashed}, {"color", style::color::green}}))
    .save(base_path + "compare.png");

    plots::PlotDataset(I_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}}))
        .plot(C_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_aa"}, {"linestyle", style::line::dashed}, {"color", style::color::green}}))
        .plot(foxs, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "foxs"}, {"color", style::color::black}}))
    .save(base_path + "vacuum.png");

    plots::PlotDataset(I_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}}))
        .plot(I_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::green}}))
        .plot(I_ww, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ww"}, {"color", style::color::brown}}))
        .plot(foxs_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aa"}, {"color", style::color::red}}))
        .plot(foxs_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xx"}, {"color", style::color::green}}))
        .plot(foxs_ww, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ww"}, {"color", style::color::brown}}))
    .save(base_path + "compare_foxs.png");

    plots::PlotDataset(I_ax, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ax"}, {"color", style::color::blue}}))
        .plot(I_aw, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aw"}, {"color", style::color::pink}}))
        .plot(I_wx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xw"}, {"color", style::color::purple}}))
        .plot(foxs_ax, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_ax"}, {"color", style::color::blue}}))
        .plot(foxs_aw, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_aw"}, {"color", style::color::pink}}))
        .plot(foxs_wx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xw"}, {"color", style::color::purple}}))
    .save(base_path + "compare_foxs_cross.png");

    plots::PlotDataset(foxs_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xx"}, {"color", style::color::green}}))
        .plot(C_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_xx"}, {"linestyle", style::line::dashed}, {"color", style::color::blue}}))
        .plot(I_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::green}}))
    .save(base_path + "exv.png");

    for (unsigned int i = 0; i < C_xx.size(); ++i) {
        C_aa.y(i) /= I_aa.interpolate_x(C_aa.x(i), 1);
        C_xx.y(i) /= I_xx.interpolate_x(C_xx.x(i), 1);
    }
    // C_xx.normalize();
    // C_aa.normalize();

    for (unsigned int i = 0; i < foxs_xx.size(); ++i) {
        foxs_aa.y(i) /= I_aa.interpolate_x(foxs_aa.x(i), 1);
        foxs_xx.y(i) /= I_xx.interpolate_x(foxs_xx.x(i), 1);
    }
    // foxs_aa.normalize();
    // foxs_xx.normalize();

    plots::PlotDataset(foxs_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"linestyle", style::line::dashed}, {"ylabel", "I"}, {"legend", "foxs_xx"}, {"color", style::color::green}}))
        .plot(C_xx, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_xx"}, {"linestyle", style::line::dashed}, {"color", style::color::blue}}))
        .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(base_path + "exv_normalized.png");

    plots::PlotDataset(foxs_aa, plots::PlotOptions({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}}))
        .plot(C_aa, plots::PlotOptions({{"xlabel", "q"}, {"linewidth", 2}, {"ylabel", "I"}, {"legend", "CRYSOL_aa"}, {"linestyle", style::line::dashed}, {"color", style::color::green}}))
        .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    .save(base_path + "vacuum_normalized.png");

    // for (unsigned int i = 0; i < foxs_xx.size(); ++i) {foxs_xx.y(i) /= foxs_aa.y(i);}
    // for (unsigned int i = 0; i < I_xx.size(); ++i) {I_xx.y(i) /= I_aa.y(i);}

    // plots::PlotDataset(foxs_xx, plots::PlotOptions({{"logy", false}, {"logx", false}, {"xlimits", std::vector<double>{0, 0.5}}, {"ylimits", std::vector<double>{0.8, 1.2}}}))
    //     .plot(I_xx)
    //     .hline(1, plots::PlotOptions(style::draw::line, {{"linestyle", style::line::dashed}, {"color", style::color::black}}))
    // .save(base_path + "xx_div_aa.png");
}