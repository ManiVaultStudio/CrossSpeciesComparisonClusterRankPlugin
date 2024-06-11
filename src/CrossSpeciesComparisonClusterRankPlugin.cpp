#include "CrossSpeciesComparisonClusterRankPlugin.h"

#include "ChartWidget.h"

#include <DatasetsMimeData.h>

#include <vector>
#include <random>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QVariantList>
#include <QVariantMap>
#include <QMimeData>
#include <QDebug>


#include "lib/JSONnlohmann/json.hpp"
#include "lib/Clustering/fastcluster.h"
#include "lib/JSONnlohmann/json.hpp"
#include "lib/Clustering/fastcluster.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <QJsonDocument>
#include <QJsonObject>
#include <QString>
#include <execution>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <stack>
#include "lib/Distance/annoylib.h"
#include "lib/Distance/kissrandom.h"

Q_PLUGIN_METADATA(IID "studio.manivault.CrossSpeciesComparisonClusterRankPlugin")

using namespace mv;

CrossSpeciesComparisonClusterRankPlugin::CrossSpeciesComparisonClusterRankPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _chartWidget(nullptr),
    _dropWidget(nullptr),
    _currentDataSet(nullptr),
    _settingsAction(*this)
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

   
    auto mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 0, 0, 0);
    mainLayout->setSpacing(0);

    auto mainOptionsLayout = new QHBoxLayout();
    mainOptionsLayout->setSpacing(0);
    mainOptionsLayout->setContentsMargins(0, 0, 0, 0);
    auto extraOptionsGroup = new VerticalGroupAction(this, "Settings");

    extraOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("cog"));
    extraOptionsGroup->addAction(&_settingsAction.getHierarchyTopClusterDataset());
    extraOptionsGroup->addAction(&_settingsAction.getHierarchyMiddleClusterDataset());
    extraOptionsGroup->addAction(&_settingsAction.getHierarchyBottomClusterDataset());
    extraOptionsGroup->addAction(&_settingsAction.getSelectedClusterNames());
    extraOptionsGroup->addAction(&_settingsAction.getOptionSelectionAction());
    extraOptionsGroup->addAction(&_settingsAction.getSpeciesNamesDataset());
    extraOptionsGroup->addAction(&_settingsAction.getMainPointsDataset());
    extraOptionsGroup->addAction(&_settingsAction.getEmbeddingDataset());
    extraOptionsGroup->addAction(&_settingsAction.getReferenceTreeDataset());
    extraOptionsGroup->addAction(&_settingsAction.getGeneNamesConnection());
    extraOptionsGroup->addAction(&_settingsAction.getFilterTreeDataset());
    
    auto mainOptionsGroup = new HorizontalGroupAction(this, "Trigger");
    mainOptionsGroup->setIcon(Application::getIconFont("FontAwesome").getIcon("play"));
    mainOptionsGroup->addAction(&_settingsAction.getCreatePointSelectTree());

    mainOptionsLayout->addWidget(mainOptionsGroup->createWidget(&getWidget()), 2);
    mainOptionsLayout->addWidget(extraOptionsGroup->createCollapsedWidget(&getWidget()), 1);

    mainLayout->addLayout(mainOptionsLayout);
    mainLayout->addWidget(_chartWidget, 1);

    getWidget().setLayout(mainLayout);
   

    const auto mainPointsDatasetUpdate = [this]() -> void
        {
            _currentDataSet = _settingsAction.getMainPointsDataset().getCurrentDataset();
            if (_currentDataSet.isValid())
            {
                auto children = _currentDataSet->getChildren();
                for (auto child : children)
                {
                    child->setGroupIndex(1);
                }
            }
        };
    connect(&_settingsAction.getMainPointsDataset(), &DatasetPickerAction::currentIndexChanged, this, mainPointsDatasetUpdate);

    const auto embeddingDatasetUpdate = [this]() -> void
        {
            _embeddingDataset = _settingsAction.getEmbeddingDataset().getCurrentDataset();
        };
    connect(&_settingsAction.getEmbeddingDataset(), &DatasetPickerAction::currentIndexChanged, this, embeddingDatasetUpdate);

    const auto clusterSelectionFromPopulationPyramidDatasetChange = [this]() -> void
        {
            qDebug() << "Item selected in Population Pyramid";
        };
    connect(&_embeddingDataset, &mv::Dataset<Points>::dataSelectionChanged, this, clusterSelectionFromPopulationPyramidDatasetChange);

    const auto getCreatePointSelectTreeUpdate = [this]() -> void
        {
            
            auto clusterDataset = _settingsAction.getHierarchyBottomClusterDataset().getCurrentDataset();
            auto pointsDataset = _settingsAction.getMainPointsDataset().getCurrentDataset();
            //auto speciesDataset= _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
            std::vector<std::uint32_t> selectedIndices= pointsDataset->getSelectionIndices();
            
            if (clusterDataset.isValid() && pointsDataset.isValid())
            {

            }

            if (_settingsAction.getFilterTreeDataset().getCurrentDataset().isValid() && _settingsAction.getSpeciesNamesDataset().getCurrentDataset().isValid() && selectedIndices.size() > 0)
            {


                auto speciesDataset = mv::data().getDataset<Clusters>(_settingsAction.getSpeciesNamesDataset().getCurrentDataset().getDatasetId());
                std::map<QString, int> speciesSelectedIndicesCounter;
                auto speciesData = speciesDataset->getClusters();
                for (auto species : speciesData)
                {
                    auto speciesNameKey = species.getName();
                    int speciescellCountValue = 0;
                    auto indices = species.getIndices();
                    speciescellCountValue = std::count_if(selectedIndices.begin(), selectedIndices.end(), [&](const auto& element) {
                        return std::find(indices.begin(), indices.end(), element) != indices.end();
                        });

                    speciesSelectedIndicesCounter.insert({ speciesNameKey, speciescellCountValue });
                }

                if (speciesSelectedIndicesCounter.size() > 0)
                {
                    // qDebug() << "CrossSpeciesComparisonClusterRankPlugin::publishSelection: Send selection to core";
                    QJsonObject valueStringReference = createJsonTree(speciesSelectedIndicesCounter);
                    //print speciesSelectedIndicesCounter
                    /*
                    qDebug() << "*******************";
                    int i = 1;
                    for (auto& [key, value] : speciesSelectedIndicesCounter)
                    {
                        qDebug() << i<< key << " : " << value;
                        i++;
                    }
                    qDebug() << "*******************";
                    */

                    auto mainTreeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilterTreeDataset().getCurrentDataset().getDatasetId());
                    mainTreeDataset->setTreeData(valueStringReference);

                    events().notifyDatasetDataChanged(mainTreeDataset);
                    _settingsAction.getGeneNamesConnection().setString("");
                }

            }
            else
            {
                qDebug() << "Datasets not valid";
            }

        };
    connect(&_settingsAction.getCreatePointSelectTree(), &TriggerAction::triggered, this, getCreatePointSelectTreeUpdate);


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

    //qDebug() << "CrossSpeciesComparisonClusterRankPlugin::loadData: Load data set from ManiVault core";

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

    //qDebug() << "CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart: Prepare payload";
    _dropWidget->setShowDropIndicator(false);
    _currentDataSet;
    QVariantList dataForChart = {
    createNode("All", "white", {
        createNode("Non-Neuronal", "#7fffff", {
            createNode("Micro-PVM", "#94af97", {
                createLeaf("Microglia/PVM", "#94af97", 20566)
            }),
            createNode("VLMC", "#697255", {
                createLeaf("No Agreement", "#C0C0C0", 133072),
                createLeaf("VLMC", "#697255", 7968),
                createLeaf("Unknown", "#C0C0C0", 50)
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
    //                            QVariantMap{{"name", "AspectRatioBanker"}, {"color", "#C0C0C0"}, {"value", 7074}},
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
       // qDebug() << "CrossSpeciesComparisonClusterRankPlugin::convertDataAndUpdateChart: Send data from Qt cpp to D3 js";
        emit _chartWidget->getCommunicationObject().qt_js_setDataAndPlotInJS(jsonString);
    }
    
}

QJsonObject CrossSpeciesComparisonClusterRankPlugin::createJsonTree(std::map<QString, int> speciesSelectedIndicesCounter)
{
    QJsonObject valueStringReference;

    {
    
        std::vector<float> numbers;
        std::vector<QString> leafnames;

        const std::unordered_map<std::string, int> clusteringTypeMap = {
    {"Complete", HCLUST_METHOD_COMPLETE},
    {"Average", HCLUST_METHOD_AVERAGE},
    {"Median", HCLUST_METHOD_MEDIAN}
        };
        std::string clusteringTypecurrentText = "Single";  //"Single","Complete", "Average","Median"
        int opt_method = clusteringTypeMap.count(clusteringTypecurrentText) ? clusteringTypeMap.at(clusteringTypecurrentText) : HCLUST_METHOD_SINGLE;

        auto numOfLeaves = speciesSelectedIndicesCounter.size();
        double* distmat = new double[(numOfLeaves * (numOfLeaves - 1)) / 2];

        //iterate std::map<QString, int> speciesSelectedIndicesCounter and populate numbers
        for (auto& [key, value] : speciesSelectedIndicesCounter)
        {
            numbers.push_back(value);
            leafnames.push_back(key);
        }





        distmat = _settingsAction.condensedDistanceMatrix(numbers);

        int* merge = new int[2 * (numOfLeaves - 1)];
        double* height = new double[numOfLeaves - 1];
        hclust_fast(numOfLeaves, distmat, opt_method, merge, height);




        std::string newick = _settingsAction.mergeToNewick(merge, numOfLeaves);
        int totalChars = newick.length();
        //append a ";" to the end of the newick string
        //qDebug()<<"Newick format: " << QString::fromStdString(newick);
        newick += ';';
        //std::cout << "Newick format: " << newick << std::endl;
        QString formattedtree = _settingsAction.createJsonTreeFromNewick(QString::fromStdString(newick),leafnames);
        valueStringReference = QJsonDocument::fromJson(formattedtree.toUtf8()).object();
    delete[] distmat;
    delete[] merge;
    delete[] height;
}
    

    return valueStringReference;
}



void CrossSpeciesComparisonClusterRankPlugin::publishSelection(const std::vector<QString>& selectedIDs)
{
    auto clusterDataset= _settingsAction.getHierarchyBottomClusterDataset().getCurrentDataset();
    auto pointsDataset= _settingsAction.getMainPointsDataset().getCurrentDataset();
    //auto speciesDataset= _settingsAction.getSpeciesNamesDataset().getCurrentDataset();
    std::vector<std::uint32_t> selectedIndices;

    if (clusterDataset.isValid() && pointsDataset.isValid())
    {
        bool isValid = false;
        isValid = clusterDataset->getParent().getDatasetId() == pointsDataset->getId();//&& speciesDataset->getParent().getDatasetId() == pointsDataset->getId();


        QList<QString> list;
        for (const auto& id : selectedIDs) {
            list.append(id);
        }

        QVariant variant = QVariant::fromValue(list);
        _settingsAction.getSelectedClusterNames().setVariant(variant);

        if (isValid)
        {
           
            if (!selectedIDs.empty())
            {

                auto rawData = mv::data().getDataset < Clusters>(clusterDataset.getDatasetId());
                auto clusters = rawData->getClusters();
                if (clusters.size() > 0)
                {
                    for (auto cluster : clusters)
                    {
                        auto clusterName = cluster.getName();
                        //if selectedIDs contains cluster.getName()
                        if (std::find(selectedIDs.begin(), selectedIDs.end(), clusterName) != selectedIDs.end())
                        {
                            auto clusterIndices = cluster.getIndices();
                            selectedIndices.insert(selectedIndices.end(), clusterIndices.begin(), clusterIndices.end());
                        }




                    }
                }

            }
            pointsDataset->setSelectionIndices(selectedIndices);
            mv::events().notifyDatasetDataSelectionChanged(pointsDataset);
        }
        else
        {
            qDebug() << "Datasets not valid";
        }

    }

    else
    {
        qDebug() << "Datasets not valid";
    }
    

    
    if (_settingsAction.getFilterTreeDataset().getCurrentDataset().isValid() && _settingsAction.getSpeciesNamesDataset().getCurrentDataset().isValid() && selectedIndices.size()>0 )
    {
        
        
        auto speciesDataset = mv::data().getDataset<Clusters>(_settingsAction.getSpeciesNamesDataset().getCurrentDataset().getDatasetId());
        std::map<QString, int> speciesSelectedIndicesCounter;
        auto speciesData = speciesDataset->getClusters();
        for (auto species : speciesData)
        {
            auto speciesNameKey= species.getName();
            int speciescellCountValue = 0;
            auto indices= species.getIndices();
            speciescellCountValue = std::count_if(selectedIndices.begin(), selectedIndices.end(), [&](const auto& element) {
                return std::find(indices.begin(), indices.end(), element) != indices.end();
                });
            
            speciesSelectedIndicesCounter.insert({ speciesNameKey, speciescellCountValue });
        }

        if ( speciesSelectedIndicesCounter.size()>0)
        {
           // qDebug() << "CrossSpeciesComparisonClusterRankPlugin::publishSelection: Send selection to core";
            QJsonObject valueStringReference = createJsonTree(speciesSelectedIndicesCounter);
            //print speciesSelectedIndicesCounter
            /*
            qDebug() << "*******************";
            int i = 1;
            for (auto& [key, value] : speciesSelectedIndicesCounter)
            {
                qDebug() << i<< key << " : " << value;
                i++;
            }
            qDebug() << "*******************";
            */

            auto mainTreeDataset = mv::data().getDataset<CrossSpeciesComparisonTree>(_settingsAction.getFilterTreeDataset().getCurrentDataset().getDatasetId());
            mainTreeDataset->setTreeData(valueStringReference);

            events().notifyDatasetDataChanged(mainTreeDataset);
            _settingsAction.getGeneNamesConnection().setString("");
        }

    }
    else
    {
        qDebug() << "Datasets not valid";
    }


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
/*
namespace  local
{
    void get_recursive_cluster_tree(QStandardItem* item, Dataset<DatasetImpl> currentDataset, const QVector<QString>& hierarchy, qsizetype h, bool firstTime, bool intersection = true, const std::vector<uint32_t>& indices = {})
    {
        if (h >= hierarchy.size())
            return;

        auto childDatasets = currentDataset->getChildren({ ClusterType });
        for (qsizetype c = 0; c < childDatasets.size(); ++c)
        {
            if (childDatasets[c]->getGuiName() == hierarchy[h])
            {
                mv::Dataset<Clusters> clusterData = childDatasets[c];
                auto clusters = clusterData->getClusters();


                if (intersection && !firstTime)
                {
                    QSet<QString> clusterNames;
                    for (auto cluster : clusters)
                    {
                        QString name = cluster.getName();
                        clusterNames.insert(name);
                    }
                    for (qsizetype row = item->rowCount() - 1; row >= 0; --row)
                    {
                        QStandardItem* child = item->child(row, 0);
                        if (!clusterNames.contains(child->text()))
                        {
                            item->removeRow(row);
                        }
                    }
                }

                for (auto cluster : clusters)
                {
                    QString name = cluster.getName();
                    QStandardItem* correspondingItem = nullptr;

                    if (!firstTime)
                    {
                        for (qsizetype row = 0; row < item->rowCount(); ++row)
                        {
                            QStandardItem* child = item->child(row, 0);
                            if (child->text() == name)
                            {
                                correspondingItem = child;
                                break;
                            }
                        }
                        if (intersection && (correspondingItem == nullptr))
                            continue;;
                    }



                    {
                        std::vector<uint32_t> clusterIndices = cluster.getIndices();
                        std::sort(clusterIndices.begin(), clusterIndices.end());
                        std::vector<uint32_t> intersectingIndices;
                        if (h == 0)
                            intersectingIndices = clusterIndices;
                        else
                        {
                            std::set_intersection(indices.cbegin(), indices.cend(), clusterIndices.cbegin(), clusterIndices.cend(), std::back_inserter(intersectingIndices));
                        }
                        if (!intersectingIndices.empty())
                        {
                            if (correspondingItem == nullptr)
                            {
                                QPixmap pixmap(16, 16);
                                pixmap.fill(cluster.getColor());
                                correspondingItem = new QStandardItem(pixmap, name);
                                correspondingItem->setData(h, Qt::UserRole);

                                item->appendRow(correspondingItem);
                            }
                            get_recursive_cluster_tree(correspondingItem, currentDataset, hierarchy, h + 1, firstTime, intersection, intersectingIndices);
                        }

                    }
                }
                break;
            }
        }
    }
}


void CrossSpeciesComparisonClusterRankPlugin::createHierarchy(qsizetype index, const Dataset<DatasetImpl>& dataset)
{
    
    const qsizetype NrOfDatasets = _selectedDatasetsAction.size();


    QVector<QString> hierarchy = { "class", "subclass","cross_species_cluster" };

    for (qsizetype i = 0; i < _selectedDatasetsAction.size(); ++i)
    {
        if (!_selectedDatasetsAction.getDataset(i).isValid())
            return;
    }
    _model.clear();
    //  _model.invisibleRootItem()->setData("Hierarchy", Qt::DisplayRole);

    QStringList finalLevelItems;
    for (qsizetype i = 0; i < _selectedDatasetsAction.size(); ++i)
    {

        Dataset<DatasetImpl> currentDataset = _selectedDatasetsAction.getDataset(i).get();


        auto childDatasets = currentDataset->getChildren({ ClusterType });
        for (qsizetype c = 0; c < childDatasets.size(); ++c)
        {
            if (childDatasets[c]->getGuiName() == hierarchy.last())
            {
                Dataset<Clusters> clusterData = childDatasets[c];
                auto clusters = clusterData->getClusters();

                for (auto cluster : clusters)
                {
                    QString name = cluster.getName();
                    finalLevelItems.append(name);
                }
            }
        }


        local::get_recursive_cluster_tree(_model.invisibleRootItem(), currentDataset, hierarchy, 0, (i == 0), true);

    }


    _selectedOptionsAction.setOptions(finalLevelItems);
    if (finalLevelItems.count())
        _selectedOptionsAction.selectOption(finalLevelItems.first());

   
}
 */



QString CrossSpeciesComparisonClusterRankPlugin::getCurrentDataSetID() const
{
    if (_currentDataSet.isValid())
        return _currentDataSet->getId();
    else
        return QString{};
}
void CrossSpeciesComparisonClusterRankPlugin::fromVariantMap(const QVariantMap& variantMap)
{
    ViewPlugin::fromVariantMap(variantMap);

    mv::util::variantMapMustContain(variantMap, "CSCCR:Cross-Species Comparison Cluster Rank Settings");
    _settingsAction.fromVariantMap(variantMap["CSCCR:Cross-Species Comparison Cluster Rank Settings"].toMap());


}

QVariantMap CrossSpeciesComparisonClusterRankPlugin::toVariantMap() const
{
    QVariantMap variantMap = ViewPlugin::toVariantMap();

    _settingsAction.insertIntoVariantMap(variantMap);

    return variantMap;
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
