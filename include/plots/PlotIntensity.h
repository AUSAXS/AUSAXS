#pragma once

class TPad;

#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>

#include <memory>
#include <string>

namespace plots {
	class PlotIntensity : public Plot {
		public:
			/**
			 * @brief Copy constructor.
			 * 
			 * @param d The ScatteringHistogram to be plotted. 
			 */
			PlotIntensity(const hist::ScatteringHistogram& d, int color = kBlack);

			/**
			 * @brief Move constructor.
			 * 
			 * @param d The ScatteringHistogram to be plotted. 
			 */
			PlotIntensity(hist::ScatteringHistogram&& d, int color = kBlack);

			/**
			 * @brief Constructor.
			 * 
			 * @param d The dataset to be plotted.
			 */
			PlotIntensity(const SimpleDataset& d);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensity() override;

			/**
			 * @brief Plot an additional data set as points. 
			 */
			void plot_intensity(const SimpleDataset& data);

			/**
			 * @brief Plot the result of a fit. 
			 */
			void plot_intensity(const std::shared_ptr<Fit> fit, const PlotOptions& options);

			void plot_guinier_approx();

			void save(std::string path) const override;
			
			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const hist::ScatteringHistogram& h, std::string path);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const SimpleDataset& d, std::string path);

		private:
			const hist::ScatteringHistogram d;
			std::shared_ptr<TCanvas> canvas;
			std::unique_ptr<TPad> linpad;
			std::unique_ptr<TPad> logpad;
			Limit limits;

			void initial_intensity_plot(int color);

			void initial_intensity_plot(const SimpleDataset& data);

			void prepare_canvas();
	};
}