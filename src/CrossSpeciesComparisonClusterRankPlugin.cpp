#include "CrossSpeciesComparisonClusterRankPlugin.h"

#include "ChartWidget.h"

#include <DatasetsMimeData.h>

#include <vector>
#include <random>
#include <ClusterData/ClusterData.h>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QVariantList>
#include <QVariantMap>
#include <QMimeData>
#include <QDebug>

Q_PLUGIN_METADATA(IID "studio.manivault.CrossSpeciesComparisonClusterRankPlugin")

using namespace mv;

CrossSpeciesComparisonClusterRankPlugin::CrossSpeciesComparisonClusterRankPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _chartWidget(nullptr),
    _dropWidget(nullptr),
    _currentDataSet(nullptr)
{
}

void CrossSpeciesComparisonClusterRankPlugin::init()
{
    getWidget().setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

    // Create layout
    auto layout = new QVBoxLayout();
    layout->setContentsMargins(0, 0, 0, 0);

    // Create chart widget and set html contents of webpage 
    _chartWidget = new ChartWidget(this);
    _chartWidget->setPage(":CrossSpeciesComparisonClusterRank_chart/icicle_chart.html", "qrc:/CrossSpeciesComparisonClusterRank_chart/");

    // Add widget to layout
    layout->addWidget(_chartWidget);

    // Apply the layout
    getWidget().setLayout(layout);

    // Instantiate new drop widget: See CrossSpeciesComparisonClusterRank for details
    _dropWidget = new DropWidget(_chartWidget);
    _dropWidget->setDropIndicatorWidget(new DropWidget::DropIndicatorWidget(&getWidget(), "No data loaded", "Drag the Point data which has child cluster datasets to this view"));

    _dropWidget->initialize([this](const QMimeData* mimeData) -> DropWidget::DropRegions {

        // A drop widget can contain zero or more drop regions
        DropWidget::DropRegions dropRegions;

        const auto datasetsMimeData = dynamic_cast<const DatasetsMimeData*>(mimeData);


        if (datasetsMimeData == nullptr)
            return dropRegions;

        if (datasetsMimeData->getDatasets().count() > 1)
            return dropRegions;

        const auto dataset = datasetsMimeData->getDatasets().first();
        
        auto children = dataset->getChildren({ ClusterType }, false);
        if (children.size() < 2)
            return dropRegions;
        
        const auto datasetGuiName = dataset->text();
        const auto datasetId = dataset->getId();
        const auto dataType = dataset->getDataType();
        const auto dataTypes = DataTypes({ PointType });


        if (dataTypes.contains(dataType) ) {

            if (datasetId == getCurrentDataSetID()) {
                dropRegions << new DropWidget::DropRegion(this, "Warning", "Data already loaded", "exclamation-circle", false);
            }
            else {
                auto candidateDataset = mv::data().getDataset<Points>(datasetId);

                dropRegions << new DropWidget::DropRegion(this, "Points", QString("Visualize %1 as cluster rank").arg(datasetGuiName), "map-marker-alt", true, [this, candidateDataset]() {
                    loadData({ candidateDataset });
                   
                    });

            }
        }
        else {
            dropRegions << new DropWidget::DropRegion(this, "Incompatible data", "This type of data is not supported", "exclamation-circle", false);
        }

        return dropRegions;
        });
    convertDataAndUpdateChart();
    // update data when data set changed
    connect(&_currentDataSet, &Dataset<Points>::dataChanged, this, &CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart);

    // Update the selection (coming from PCP) in core
    connect(&_chartWidget->getCommunicationObject(), &ChartCommObject::passSelectionToCore, this, &CrossSpeciesComparisonClusterRankPlugin::publishSelection);
}

void CrossSpeciesComparisonClusterRankPlugin::loadData(const mv::Datasets& datasets)
{
    // Exit if there is nothing to load
    if (datasets.isEmpty())
        return;

    qDebug() << "CrossSpeciesComparisonClusterRankPlugin::loadData: Load data set from ManiVault core";

    // Load the first dataset, changes to _currentDataSet are connected with convertDataAndUpdateChart
    _currentDataSet = datasets.first();
    events().notifyDatasetDataChanged(_currentDataSet);
}
QVariantMap createNode(const QString& name, const QString& color, const QVariantList& children = QVariantList())
{
    QVariantMap node;
    node.insert("name", name);
    node.insert("color", color);
    if (!children.isEmpty())
        node.insert("children", children);
    return node;
}

QVariantMap createLeaf(const QString& name, const QString& color, int value)
{
    QVariantMap leaf;
    leaf.insert("name", name);
    leaf.insert("color", color);
    leaf.insert("value", value);
    return leaf;
}

void CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart()
{
    //if (!_currentDataSet.isValid())
       // return;

    qDebug() << "CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart: Prepare payload";
    _dropWidget->setShowDropIndicator(false);
    _currentDataSet;
    QVariantList dataForChart = {
    createNode("All", "white", {
        createNode("Non-Neuronal", "#7fffff", {
            createNode("Micro-PVM", "#94af97", {
                createLeaf("Microglia/PVM", "#94af97", 20566)
            }),
            createNode("VLMC", "#697255", {
                createLeaf("No Agreement", "silver", 133072),
                createLeaf("VLMC", "#697255", 7968),
                createLeaf("Unknown", "silver", 50)
            }),
            createNode("Endo", "#8d6c62", {
                createLeaf("Endo", "#8d6c62", 9123)
            }),
            createNode("Oligo", "#53776c", {
                createLeaf("Oligo_1", "#88a19a", 63386),
                createLeaf("Oligo_2", "#3a544c", 90484)
            }),
            createNode("Astro", "#665c47", {
                createLeaf("Astro_1", "#958f80", 5583),
                createLeaf("Astro_2", "#484132", 60795)
            }),
            createNode("OPC", "#374a45", {
                createLeaf("OPC", "#374a45", 28445)
            }),
        }),
        createNode("GABAergic", "#ff7f7f", {
            createNode("Sst Chodl", "#b1b10c", {
                createLeaf("Sst Chodl", "#b1b10c", 5133)
            }),
            createNode("Sst", "#ff9900", {
                createLeaf("Sst_1", "#ffb84f", 11487),
                createLeaf("Sst_2", "#ffb342", 12334),
                createLeaf("Sst_3", "#ffb03c", 74776),
                createLeaf("Sst_4", "#ffae36", 15488),
                createLeaf("Sst_5", "#ffac30", 6547),
                createLeaf("Sst_6", "#ffa92a", 9323),
                createLeaf("Sst_7", "#ffa724", 1607)
            }),
            createNode("Pvalb", "#d93137", {
                createLeaf("Pvalb_1", "#e47175", 143140),
                createLeaf("Pvalb_2", "#e2686c", 10764)
            }),
            createNode("Chandelier", "#f641a8", {
                createLeaf("Chandelier", "#f641a8", 16199)
            }),
            createNode("Vip", "#a45fbf", {
                createLeaf("Vip_1", "#c091d3", 15946),
                createLeaf("Vip_2", "#bd8cd1", 42074),
                createLeaf("Vip_3", "#ba87cf", 26303),
                createLeaf("Vip_4", "#b883cd", 10681)

            }),
            createNode("Sncg", "#df70ff", {
                createLeaf("Sncg_1", "#e99cff", 7342),
                createLeaf("Sncg_2", "#e691ff", 4291),
                createLeaf("Sncg_3", "#e386ff", 4459),
                createLeaf("Sncg_4", "#e17bff", 13487),
            }),
            createNode("Lamp5", "#da808c", {
                createLeaf("Lamp5_1", "#e5a7b0", 4557),
                createLeaf("Lamp5_2", "#e097a1", 40821),
                createLeaf("Lamp5_3", "#dc8793", 22716),
                createLeaf("Lamp5_4", "#cd7883", 6217),
                createLeaf("Lamp5_5", "#b36972", 16354)
            })
        }),
        createNode("Glutamatergic", "#bfff7f", {
            createNode("L6 IT", "#a19922", {
                createLeaf("L6 IT_1", "#beb867", 24546),
                createLeaf("L6 IT_2", "#716c18", 136534),
                createLeaf("L6 IT_3", "#ffd700", 17004)
            }),
            createNode("L2/3 IT", "#b1ec30", {
                createLeaf("L2/3 IT", "#b1ec30", 486221)
            }),
            createNode("L5 ET", "#0d5b78", {
                createLeaf("L5 ET_1", "#588ea2", 52299),
                createLeaf("L5 ET_2", "#0d5b78", 313)

            }),
            createNode("L5 IT", "#50b2ad", {
                createLeaf("L5 IT_1", "#86cac6", 316754),
                createLeaf("L5 IT_2", "#67bcb7", 191918),
                createLeaf("L5 IT_3", "#beb867", 18024)
            }),
            createNode("L6b", "#7044aa", {
                createLeaf("L6b", "#7044aa", 61326)
            }),
            createNode("L6 CT", "#2d8cb8", {
                createLeaf("L6 CT_1", "#6eb0ce", 77283),
                createLeaf("L6 CT_2", "#4d9ec3", 119532)
            }),
            createNode("L6 IT Car3", "#5100ff", {
                createLeaf("MC", "#5100ff", 7074)
            }),
            createNode("L5/6 NP", "#3e9e64", {
                createLeaf("L5/6 NP", "#3e9e64", 47918)
            })
        })
    })
    };
    //QVariantList dataForChart = {
    //    QVariantMap{
    //        {"name", "flare"},
    //        {"color", "pink"},
    //        {"children", QVariantList{
    //            QVariantMap{
    //                {"name", "analytics"},
    //                {"color", "blue"},
    //                {"children", QVariantList{
    //                    QVariantMap{
    //                        {"name", "cluster"},
    //                        {"color", "green"},
    //                        {"children", QVariantList{
    //                            QVariantMap{{"name", "AgglomerativeCluster"}, {"color", "red"}, {"value", 3938}},
    //                            QVariantMap{{"name", "CommunityStructure"}, {"color", "orange"}, {"value", 3812}},
    //                            QVariantMap{{"name", "HierarchicalCluster"}, {"color", "yellow"}, {"value", 6714}},
    //                            QVariantMap{{"name", "MergeEdge"}, {"color", "purple"}, {"value", 743}},
    //                        }},
    //                    },
    //                    QVariantMap{
    //                        {"name", "graph"},
    //                        {"color", "cyan"},
    //                        {"children", QVariantList{
    //                            QVariantMap{{"name", "BetweennessCentrality"}, {"color", "magenta"}, {"value", 3534}},
    //                            QVariantMap{{"name", "LinkDistance"}, {"color", "lime"}, {"value", 5731}},
    //                            QVariantMap{{"name", "MaxFlowMinCut"}, {"color", "teal"}, {"value", 7840}},
    //                            QVariantMap{{"name", "ShortestPaths"}, {"color", "maroon"}, {"value", 5914}},
    //                            QVariantMap{{"name", "SpanningTree"}, {"color", "navy"}, {"value", 3416}},
    //                        }},
    //                    },
    //                    QVariantMap{
    //                        {"name", "optimization"},
    //                        {"color", "olive"},
    //                        {"children", QVariantList{
    //                            QVariantMap{{"name", "AspectRatioBanker"}, {"color", "silver"}, {"value", 7074}},
    //                        }},
    //                    },
    //                }},
    //            },
    //        }},
    //    },
    //};

    QJsonDocument doc = QJsonDocument::fromVariant(QVariant::fromValue(dataForChart));
    QString jsonString = doc.toJson(QJsonDocument::Compact);
    
    if (jsonString!="")
    {
        qDebug() << "CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart: Send data from Qt cpp to D3 js";
        emit _chartWidget->getCommunicationObject().qt_js_setDataAndPlotInJS(jsonString);
    }
    
}

void CrossSpeciesComparisonClusterRankPlugin::publishSelection(const std::vector<QString>& selectedIDs)
{
    
    qDebug() << "\nSelectedIDs: ";
    for (const auto& id : selectedIDs)
    {
        qDebug() << id+" ,";
    }
    qDebug() << "\n";
    //// ask core for the selection set for the current data set
    //auto selectionSet = _currentDataSet->getSelection<Points>();
    //auto& selectionIndices = selectionSet->indices;

    //// clear the selection and add the new points
    //selectionIndices.clear();
    //selectionIndices.reserve(_currentDataSet->getNumPoints());
    //for (const auto id : selectedIDs) {
    //    selectionIndices.push_back(id);
    //}

    //// notify core about the selection change
    //if (_currentDataSet->isDerivedData())
    //    events().notifyDatasetDataSelectionChanged(_currentDataSet->getSourceDataset<DatasetImpl>());
    //else
    //    events().notifyDatasetDataSelectionChanged(_currentDataSet);
}

QString CrossSpeciesComparisonClusterRankPlugin::getCurrentDataSetID() const
{
    if (_currentDataSet.isValid())
        return _currentDataSet->getId();
    else
        return QString{};
}


// =============================================================================
// Plugin Factory 
// =============================================================================

QIcon CrossSpeciesComparisonClusterRankPluginFactory::getIcon(const QColor& color /*= Qt::black*/) const
{
    return mv::Application::getIconFont("FontAwesome").getIcon("sitemap", color);
}

ViewPlugin* CrossSpeciesComparisonClusterRankPluginFactory::produce()
{
    return new CrossSpeciesComparisonClusterRankPlugin(this);
}

mv::DataTypes CrossSpeciesComparisonClusterRankPluginFactory::supportedDataTypes() const
{
    // This CrossSpeciesComparisonClusterRank analysis plugin is compatible with points datasets
    DataTypes supportedTypes;
    supportedTypes.append(PointType);
    return supportedTypes;
}

mv::gui::PluginTriggerActions CrossSpeciesComparisonClusterRankPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    PluginTriggerActions pluginTriggerActions;
    bool isValidDataset = false;
    const auto getPluginInstance = [this]() -> CrossSpeciesComparisonClusterRankPlugin* {
        return dynamic_cast<CrossSpeciesComparisonClusterRankPlugin*>(plugins().requestViewPlugin(getKind()));
    };

    const auto numberOfDatasets = datasets.count();

    if (numberOfDatasets == 1 && PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {

        auto dataset= datasets.first();
        auto children = dataset->getChildren({ClusterType},false);
        if (children.size()>1)
        {
            isValidDataset = true;
        }

        if (isValidDataset)
        {
            auto pluginTriggerAction = new PluginTriggerAction(const_cast<CrossSpeciesComparisonClusterRankPluginFactory*>(this), this, "Cross-Species Comparison Cluster Rank View", "Cross-Species Comparison Cluster Rank visualization", getIcon(), [this, getPluginInstance, datasets](PluginTriggerAction& pluginTriggerAction) -> void {
                for (auto dataset : datasets)
                    getPluginInstance()->loadData(Datasets({ dataset }));

                });

            pluginTriggerActions << pluginTriggerAction;
        }



    }

    return pluginTriggerActions;
}
