'''
Classes using the data to visualize the results of the analysis.

author: U.B.
'''

import pandas as pd
from dataclasses import dataclass
from pathlib import Path
from crispr_cas_evaluation.analysis.Analyzer import AnalyzerConfig, CRISPRAnalyzer
from crispr_cas_evaluation.predictions.RNARecordsAssembler import RNAmotiFoldRecordsAssembler, RNAmotiCesRecordsAssembler
from crispr_cas_evaluation.predictions.Visualization import HeatMap

#TO-DO: Add PlotConfig to CRISPRHeatmapVisualizer
@dataclass
class PlotConfig:
    '''Class for the settings of the plots.'''
    rna_type: str
    extra_info: str
    mfe_in_path: str
    folder: str
    base_plot_path: str = "./crispr_cas_evaluation/plots"

    @property
    def base_path(self) -> str:
        '''Gets the base path from the RNA type.'''
        return self.rna_type.replace(" ", "").lower()

    def shape_plot_path(self, plot_type: str) -> str:
        '''Creates the path for the shape abstraction files.'''
        return f"{self.base_plot_path}/shape_abstraction/{self.folder}/{self.base_path}_{self.mfe_in_path}_{plot_type}.jpg"

    def positional_plot_path(self, plot_type: str) -> str:
        '''Creates the path for the positional abstraction files.'''
        return f"{self.base_plot_path}/positional_abstraction/{self.folder}/{self.base_path}_{self.mfe_in_path}_{plot_type}.jpg"
    
class MotifoldPlotter:
    '''Class plotting the RNAmotiFOld results.'''
    def __init__(self, assembly: RNAmotiFoldRecordsAssembler, config: PlotConfig) -> None:
        '''Initializes a plotter object for the RNAmotiFold assemblies.'''
        self.assembly = assembly
        self.config = config

    def barchart(self) -> None:
        '''Creates and saves a barchart of the RNAmotiFold assembly.'''
        self.assembly.visualize_as_barchart(
            f"{self.config.rna_type}: Motif occurrences per sequence\n{self.config.extra_info}",
            self.config.shape_plot_path("barchart")
        )

    def violinplot(self) -> None:
        '''Creates and saves a violinplot of the RNAmotiFold assembly.'''
        self.assembly.visualize_as_violinplot(
            f"{self.config.rna_type}: Distance of mfe to lowest mfe\n{self.config.extra_info}",
            self.config.shape_plot_path("violinplot")
        )

class MoticesPlotter:
    '''Class plotting the RNAmotiFOld results.'''
    def __init__(self, assembly: RNAmotiCesRecordsAssembler, config: PlotConfig) -> None:
        '''Initializes a plotter object for the RNAHeliCes/RNAmotiCes assemblies.'''
        self.assembly = assembly
        self.config = config

    def histogram(self) -> None:
        '''Creates and saves a histogram of the RNAHeliCes/RNAmotiCes assembly.'''
        self.assembly.visualize_as_histogram(
            f"{self.config.rna_type}: Common sequences in hairpins\n{self.config.extra_info}",
            self.config.positional_plot_path("histogram")
        )

class CRISPRHeatmapVisualizer:
    '''Class visualizing the CRISPR RNA median distance of mfe for all subtypes.'''
    def __init__(self, analyzer: CRISPRAnalyzer, rna_type: str) -> None:
        '''Initializes a CRISPRHeatmapVisualizer object.'''
        self.analyzer = analyzer
        self.rna_type = rna_type

    def generate_heatmap(self, data_dict: dict, mfe_in_path: str, extra_info: str) -> None:
        '''Generates the heatmap and saves the plot in the specified path.'''
        df = pd.DataFrame(data_dict)
        heatmap = HeatMap(
            f"{self.rna_type}: Median distance of mfe\n{extra_info}",
            df
        )
        heatmap.save_plot(
            f"./crispr_cas_evaluation/plots/shape_abstraction/heatmaps/{self.rna_type.replace(" ", "").lower()}_{mfe_in_path}_heatmap.jpg"
        )

    def visualize(self) -> None:
        '''Visualizes the data of the RNAmotiFold assemblies as heatmaps.'''
        self._visualize_motifold_data()
        self._visualize_motifold_data_mfe_range()

    def _visualize_motifold_data(self) -> None:
        '''Generates heatmap for all RNAs and subtypes without mfe filtering.'''
        all_medians = {"All RNAs": self.analyzer.rna_motifold_assembly.median_distance_to_lowest_all_motifs}
        subtype_medians = {
            subtype: self.analyzer.filter_motifold_by_subtype(subtype).median_distance_to_lowest_all_motifs
            for subtype in sorted(self.analyzer.rna_motifold_assembly.unique_subtypes)
        }
        self.generate_heatmap({**all_medians, **subtype_medians}, "mfe_below0", "mfe<0 kcal/mol")

    def _visualize_motifold_data_mfe_range(self) -> None:
        '''Generates heatmap for RNAs and subtypes with MFE range filtering.'''
        motifold_mfe = self.analyzer.filter_motifold_mfe_range()
        all_medians = {"All RNAs": motifold_mfe.median_distance_to_lowest_all_motifs}
        subtype_medians = {
            subtype: self.analyzer.filter_motifold_by_subtype(subtype, mfe_range=True).median_distance_to_lowest_all_motifs
            for subtype in sorted(self.analyzer.rna_motifold_assembly.unique_subtypes)
        }
        self.generate_heatmap({**all_medians, **subtype_medians}, "mfe_range", "mfe range")

class CRISPRRNAPredictionVisualizer:
    '''Class visualizing the CRISPR RNA predictions.'''
    def __init__(self, config: AnalyzerConfig, rna_type: str) -> None:
        '''Initializes a CRISPRRNAPredictionVisualizer object.'''
        self.analyzer = CRISPRAnalyzer(config)
        self.rna_type = rna_type

    def _visualize_mfe_below_zero(self) -> None:
        '''Visualize the data for all RNAs for mfe<0'''
        config = PlotConfig(self.rna_type, "mfe<0 kcal/mol", "mfe_below0", "all_shapes")
        MotifoldPlotter(self.analyzer.rna_motifold_assembly, config).barchart()
        MotifoldPlotter(self.analyzer.rna_motifold_assembly, config).violinplot()
        MoticesPlotter(self.analyzer.rna_motices_assembly, config).histogram()

    def _visualize_mfe_range(self) -> None:
        '''Visualize the data for all RNAs within the mfe range.'''
        motifold = self.analyzer.filter_motifold_mfe_range()
        motices = self.analyzer.filter_motices_mfe_range()

        mfe_info = f"mfe range: {abs(motifold.lowest_mfe_value * 0.1):.2f} kcal/mol"
        config = PlotConfig(self.rna_type, mfe_info, "mfe_range", "all_shapes")
        MotifoldPlotter(motifold, config).barchart()
        MotifoldPlotter(motifold, config).violinplot()

        mfe_info = f"mfe range: {abs(motices.lowest_mfe_value * 0.1):.2f} kcal/mol"
        config = PlotConfig(self.rna_type, mfe_info, "mfe_range", "all_shapes")
        MoticesPlotter(motices, config).histogram()

    def visualize_all_data(self) -> None:
        '''Visualize the data for all RNAs for mfe<0 and within the mfe range.'''
        self._visualize_mfe_below_zero()
        self._visualize_mfe_range()

    def _visualize_motifold_subtypes(self) -> None:
        '''Visualizes the data of the RNAmotiFold predictions for every subtype for mfe<0 and within the mfe range.'''
        for subtype in self.analyzer.rna_motifold_assembly.unique_subtypes:
            motifold = self.analyzer.filter_motifold_by_subtype(subtype)
            config = PlotConfig(self.rna_type, f"mfe<0 kcal/mol, subtype: {subtype}", "mfe_below0", f"subtypes/{subtype}")
            MotifoldPlotter(motifold, config).barchart()
            MotifoldPlotter(motifold, config).violinplot()

            motifold_mfe = self.analyzer.filter_motifold_by_subtype(subtype, mfe_range=True)
            mfe_info = f"mfe range: {abs(motifold_mfe.lowest_mfe_value * 0.1):.2f} kcal/mol, subtype: {subtype}"
            config = PlotConfig(self.rna_type, mfe_info, "mfe_range", f"subtypes/{subtype}")
            MotifoldPlotter(motifold_mfe, config).barchart()
            MotifoldPlotter(motifold_mfe, config).violinplot()

    def _visualize_motices_subtypes(self) -> None:
        '''Visualizes the data of the RNAHeliCes/RNAmotiCes predictions for every subtype for mfe<0 and within the mfe range.'''
        for subtype in self.analyzer.rna_motices_assembly.unique_subtypes:
            motices = self.analyzer.filter_motices_by_subtype(subtype)
            config = PlotConfig(self.rna_type, "mfe<0 kcal/mol", "mfe_below0", f"subtypes/{subtype}")
            MoticesPlotter(motices, config).histogram()

            motices_mfe = self.analyzer.filter_motices_by_subtype(subtype, mfe_range=True)
            mfe_info = f"mfe range: {abs(motices_mfe.lowest_mfe_value * 0.1):.2f} kcal/mol, subtype: {subtype}"
            config = PlotConfig(self.rna_type, mfe_info, "mfe_range", f"subtypes/{subtype}")
            MoticesPlotter(motices_mfe, config).histogram()

    def visualize_subtypes(self) -> None:
        '''Visualizes the data for every subtype.'''
        self._visualize_motifold_subtypes()
        self._visualize_motices_subtypes()


    def visualize_heatmaps(self) -> None:
        '''Visualizes the heatmaps for all data and all subtypes.'''
        CRISPRHeatmapVisualizer(self.analyzer, self.rna_type).visualize()