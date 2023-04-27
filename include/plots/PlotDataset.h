#pragma once

#include <plots/Plot.h>
#include <dataset/Dataset.h>
#include <dataset/Multiset.h>

#include <concepts>

namespace plots {
	template<typename C>
	concept DatasetType = std::is_base_of_v<Dataset, C> || std::is_base_of_v<Multiset, C>;

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
			void plot(const T& data);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			template<DatasetType T>
			static void quick_plot(const T& data, const io::File& path);
		};
}