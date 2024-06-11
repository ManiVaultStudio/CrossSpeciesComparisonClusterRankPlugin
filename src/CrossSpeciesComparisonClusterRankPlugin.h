#pragma once

#include <ViewPlugin.h>

#include <Dataset.h>
#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include<CrossSpeciesComparisonTreeData.h>
#include <widgets/DropWidget.h>
#include "SettingsAction.h"
#include <QWidget>
#include <actions/HorizontalToolbarAction.h>

/** All plugin related classes are in the ManiVault plugin namespace */
using namespace mv::plugin;

/** Drop widget used in this plugin is located in the ManiVault gui namespace */
using namespace mv::gui;

/** Dataset reference used in this plugin is located in the ManiVault util namespace */
using namespace mv::util;

class ChartWidget;

/**
 * CrossSpeciesComparisonClusterRank view JS plugin class
 * 
 * This plugin showcases how a JavaScript-based visualization can be included in ManiVault 
 * Here, we use a D3 library, but other libraries like Vega-Lite follow the same pattern
 * 
 * This project:
 *  - Sets up a WebWidget, which displays an HTML webpage
 *  - Connects selections made a D3 plot with ManiVault
 * 
 * This projects does not implement selections from ManiVault to the D3 plot,
 * but such implementation follows the same form as the data-values communication
 * between cpp and JavaScript that is used here.
 *
 * @authors J. Thijssen & T. Kroes & A. Vieth
 */
class CrossSpeciesComparisonClusterRankPlugin : public ViewPlugin
{
    Q_OBJECT

public:
    /**
     * Constructor
     * @param factory Pointer to the plugin factory
     */
    CrossSpeciesComparisonClusterRankPlugin(const PluginFactory* factory);

    /** Destructor */
    ~CrossSpeciesComparisonClusterRankPlugin() override = default;
    
    /** This function is called by the core after the view plugin has been created */
    void init() override;

    /** Store a private reference to the data set that should be displayed */
    void loadData(const mv::Datasets& datasets) override;
    //void createHierarchy(qsizetype index, const Dataset<DatasetImpl>& dataset);
public slots:
    /** Converts ManiVault's point data to a json-like data structure that Qt can pass to the JS code */
    void convertDataAndUpdateChart();


private:
    /** Published selections received from the JS side to ManiVault's core */
    void publishSelection(const std::vector<QString>& selectedIDs);
    QJsonObject createJsonTree(std::map<QString, int> speciesSelectedIndicesCounter);
    QString getCurrentDataSetID() const;
    SettingsAction& getSettingsAction() { return _settingsAction; }
public: // Serialization

    /**
     * Load plugin from variant map
     * @param Variant map representation of the plugin
     */
    Q_INVOKABLE void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save plugin to variant map
     * @return Variant map representation of the plugin
     */
    Q_INVOKABLE QVariantMap toVariantMap() const override;
protected:
    ChartWidget*            _chartWidget;       // WebWidget that sets up the HTML page
    SettingsAction      _settingsAction;    // Settings action for the plugin
    DropWidget*             _dropWidget;        // Widget for drag and drop behavior
    mv::Dataset<Points>   _currentDataSet;    // Reference to currently shown data set
    //mv::Dataset<Clusters> _clusterDataset;    // Reference to the cluster dataset
    mv::Dataset<Points>   _embeddingDataset; // Reference to the low-dimensional t-SNE dataset
    //Dataset<CrossSpeciesComparisonTree>    _mainTreeDataset; // Reference to the main tree dataset
};

/**
 * CrossSpeciesComparisonClusterRank view plugin factory class
 *
 * Note: Factory does not need to be altered (merely responsible for generating new plugins when requested)
 */
class CrossSpeciesComparisonClusterRankPluginFactory : public ViewPluginFactory
{
    Q_INTERFACES(mv::plugin::ViewPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "studio.manivault.CrossSpeciesComparisonClusterRankPlugin"
                      FILE  "CrossSpeciesComparisonClusterRankPlugin.json")

public:

    /** Default constructor */
    CrossSpeciesComparisonClusterRankPluginFactory() {}

    /** Destructor */
    ~CrossSpeciesComparisonClusterRankPluginFactory() override {}
    
    /** Get plugin icon */
    QIcon getIcon(const QColor& color = Qt::black) const override;

    /** Creates an instance of the CrossSpeciesComparisonClusterRank view plugin */
    ViewPlugin* produce() override;

    /** Returns the data types that are supported by the CrossSpeciesComparisonClusterRank view plugin */
    mv::DataTypes supportedDataTypes() const override;

    /**
     * Get plugin trigger actions given \p datasets
     * @param datasets Vector of input datasets
     * @return Vector of plugin trigger actions
     */
    PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;
};
