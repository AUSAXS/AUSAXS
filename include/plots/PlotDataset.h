#pragma once

#include <plots/Plot.h>

#include <concepts>

class Dataset;
class Multiset;
namespace plots {
	template<typename C>
	concept DatasetType = std::is_base_of_v<Dataset, C> || std::is_base_of_v<Multiset, C>;

	class PlotOptions;

	/**
	 * @brief Plot the contents of a specific \class Dataset object.
	 * 		  Multiple datasets can be plotted by chaining calls to the plot() method.
	 */
	class PlotDataset : public Plot {
		public:
			/**
			 * @brief Default constructor.
			 */
			PlotDataset() {}

			/**
			 * @brief Constructor.
			 */
			template<DatasetType T>
			PlotDataset(const T& d);

			/**
			 * @brief Constructor.
			 */
			PlotDataset(const Multiset& d);

			/**
			 * @brief Destructor.
			 */
			~PlotDataset() override;

			/**
			 * @brief Plot an additional Dataset. 
			 */
			template<DatasetType T> 
			PlotDataset& plot(const T& data);

			/**
			 * @brief Plot a vertical line at the specified x coordinate.
			 */
			PlotDataset& vline(double x, const PlotOptions& options);

			/**
			 * @brief Plot a horizontal line at the specified y coordinate.
			 */
			PlotDataset& hline(double y, const PlotOptions& options);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			template<DatasetType T>
			static void quick_plot(const T& data, const io::File& path);
		};
}