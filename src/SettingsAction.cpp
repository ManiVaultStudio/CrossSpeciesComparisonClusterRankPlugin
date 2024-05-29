#include "SettingsAction.h"
#include "CrossSpeciesComparisonClusterRankPlugin.h"
#include <vector>
#include <numeric>
#include <execution>
#include <iostream>
#include <algorithm>
#include <vector>

void printMap(const std::map<QString, std::map<QString, float>>& map) {
    for (const auto& outerPair : map) {
        //std::cout << "Cluster Name: " << outerPair.first.toStdString() << std::endl;
        for (const auto& innerPair : outerPair.second) {
            //std::cout << "    Gene Name: " << innerPair.first.toStdString() << ", Expression Value: " << innerPair.second << std::endl;
        }
    }
}
float calculateMean(const std::vector<float>& v) {
    if (v.empty())
        return 0.0f;

    float sum = std::reduce(std::execution::par, v.begin(), v.end());
    float mean = sum / v.size();

    return mean;
}
QVariant createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map) {
    
    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }

    QStandardItemModel* model = new QStandardItemModel();
    model->setHorizontalHeaderItem(0, new QStandardItem("Gene"));
    model->setHorizontalHeaderItem(1, new QStandardItem("Variance"));

    int numOfSpecies = map.size();
    std::map<QString, std::map<QString, float>>::const_iterator it = map.begin();
    for (int i = 0 + 2; i < numOfSpecies + 2; i++, it++) {
        model->setHorizontalHeaderItem(i, new QStandardItem(QString("Mean_") + it->first));
    }

    for (const auto& gene : returnGeneList) {
        QList<QStandardItem*> rowItems;
        rowItems.append(new QStandardItem(gene));

        // Gather all expression values for this gene
        std::vector<float> expressionValues;
        for (const auto& species : map.at(gene)) {
            expressionValues.push_back(species.second);
        }

        // Calculate mean
        float mean = calculateMean(expressionValues);

        // Calculate variance
        float variance = 0.0f;
        for (const auto& value : expressionValues) {
            variance += std::pow(value - mean, 2);
        }
        variance /= expressionValues.size();

        rowItems.append(new QStandardItem(QString::number(variance)));

        // Add mean values for each species
        for (const auto& species : map.at(gene)) {
            rowItems.append(new QStandardItem(QString::number(species.second)));
        }

        model->appendRow(rowItems);
    }

    qRegisterMetaType<QStandardItemModel>("QStandardItemModel");
    QVariant variant = QVariant::fromValue(model);
    delete model;
    return variant;



}

QVariant  findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n) {
    
    if (map.empty() || n <= 0) {
        return QVariant();
    }

    QSet<QString> geneList;
    QStringList returnGeneList;

    for (const auto& outerPair : map) {
        //std::cout << "Species Name: " << outerPair.first.toStdString() << std::endl;

        // Convert map to vector of pairs
        std::vector<std::pair<QString, float>> geneExpressionVec(outerPair.second.begin(), outerPair.second.end());

        // Sort the vector in descending order based on the expression value
        std::sort(geneExpressionVec.begin(), geneExpressionVec.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
            });

        // Print the top 10 genes
        for (int i = 0; i < std::min(n, static_cast<int>(geneExpressionVec.size())); ++i) {
            geneList.insert(geneExpressionVec[i].first);
            //std::cout << "    Gene Name: " << geneExpressionVec[i].first.toStdString() << ", Expression Value: " << geneExpressionVec[i].second << std::endl;
        }
    }
    //convert     QSet<QString> geneList; to QStringList returnGeneList;
    for (const auto& gene : geneList) {
        returnGeneList.push_back(gene);
    }
    return createModelFromData(returnGeneList, map);
}


SettingsAction::SettingsAction(CrossSpeciesComparisonClusterRankPlugin& CrossSpeciesComparisonClusterRankPlugin) :
    WidgetAction(&CrossSpeciesComparisonClusterRankPlugin, "CrossSpeciesComparisonClusterRankPlugin Settings"),
    _crossSpeciesComparisonClusterRankPlugin(CrossSpeciesComparisonClusterRankPlugin),
    _mainPointsDataset(this, "Main Points Dataset"),
    _hierarchyTopClusterDataset(this, "Hierarchy Top Cluster Dataset"),
    _hierarchyMiddleClusterDataset(this, "Hierarchy Middle Cluster Dataset"),
    _hierarchyBottomClusterDataset(this, "Hierarchy Bottom Cluster Dataset"),
    _selectedClusterNamesVariant(this, "Selected Cluster Names Variant"),
    _filteredGeneNamesVariant(this, "Filtered Gene Names Variant"),
    _updateButtonForGeneFiltering(this, "Filter Genes Trigger"),
    _speciesNamesDataset(this, "Species Names Dataset"),
    _topNGenesFilter(this, "Top N Genes Filter"),
    _optionSelectionAction(*this)
{
    setSerializationName("CSCCR:Cross-Species Comparison Cluster Rank Settings");
    _mainPointsDataset.setSerializationName("CSCCR:MainPointsDataset");
    _hierarchyTopClusterDataset.setSerializationName("CSCCR:HierarchyTopClusterDataset");
    _hierarchyMiddleClusterDataset.setSerializationName("CSCCR:HierarchyMiddleClusterDataset");
    _hierarchyBottomClusterDataset.setSerializationName("CSCCR:HierarchyBottomClusterDataset");
    _selectedClusterNamesVariant.setSerializationName("CSCCR:SelectedClusterNamesVariant");
    _filteredGeneNamesVariant.setSerializationName("CSCCR:FilteredGeneNamesVariant");
    //_updateButtonForGeneFiltering.setSerializationName("CSCCR:UpdateButtonForGeneFiltering");
    _speciesNamesDataset.setSerializationName("CSCCR:SpeciesNamesDataset");
    _topNGenesFilter.setSerializationName("CSCCR:TopNGenesFilter");

    setText("Cross-Species Comparison Cluster Rank Settings");
    _mainPointsDataset.setToolTip("Main Points Dataset");
    _hierarchyTopClusterDataset.setToolTip("Hierarchy Top Cluster Dataset");
    _hierarchyMiddleClusterDataset.setToolTip("Hierarchy Middle Cluster Dataset");
    _hierarchyBottomClusterDataset.setToolTip("Hierarchy Bottom Cluster Dataset");
    _selectedClusterNamesVariant.setToolTip("Selected Cluster Names Variant");
    _filteredGeneNamesVariant.setToolTip("Filtered Gene Names Variant");
    _updateButtonForGeneFiltering.setToolTip("Filter Genes Trigger");
    _speciesNamesDataset.setToolTip("Species Names Dataset");
    _topNGenesFilter.setToolTip("Top N Genes Filter");
    _topNGenesFilter.initialize(1, 100, 10);


    _mainPointsDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == PointType;
        });

    _hierarchyTopClusterDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _hierarchyMiddleClusterDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });

    _hierarchyBottomClusterDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });
    _speciesNamesDataset.setFilterFunction([this](mv::Dataset<DatasetImpl> dataset) -> bool {
        return dataset->getDataType() == ClusterType;
        });

    const auto hierarchyTopClusterDatasetUpdate = [this]() -> void
        {

        };
    connect(&_hierarchyTopClusterDataset, &DatasetPickerAction::currentIndexChanged, this, hierarchyTopClusterDatasetUpdate);

    const auto hierarchyMiddleClusterDatasetUpdate = [this]() -> void
        {

        };
    connect(&_hierarchyMiddleClusterDataset, &DatasetPickerAction::currentIndexChanged, this, hierarchyMiddleClusterDatasetUpdate);



    const auto hierarchyBottomClusterDatasetUpdate = [this]() -> void
        {

        };
    connect(&_hierarchyBottomClusterDataset, &DatasetPickerAction::currentIndexChanged, this, hierarchyBottomClusterDatasetUpdate);

    const auto mainPointsDatasetUpdate = [this]() -> void
        {

        };
    connect(&_mainPointsDataset, &DatasetPickerAction::currentIndexChanged, this, mainPointsDatasetUpdate);


    const auto filteredGeneNamesVariantUpdate = [this]() -> void
        {

        };
    connect(&_filteredGeneNamesVariant, &VariantAction::changed, this, filteredGeneNamesVariantUpdate);

    const auto selectedClusterNamesVariantUpdate = [this]() -> void
        {

        };
    connect(&_selectedClusterNamesVariant, &VariantAction::changed, this, selectedClusterNamesVariantUpdate);

    const auto updateButtonForGeneFilteringUpdate = [this]() -> void
        {
            auto clusterDataset = _hierarchyBottomClusterDataset.getCurrentDataset();
            auto pointsDataset = _mainPointsDataset.getCurrentDataset();
            auto speciesDataset = _speciesNamesDataset.getCurrentDataset();
            bool isValid = false;
            _clusterNameToGeneNameToExpressionValue.clear();
            if (clusterDataset.isValid() && pointsDataset.isValid() && speciesDataset.isValid())
            {
                isValid = clusterDataset->getParent().getDatasetId() == pointsDataset->getId() && speciesDataset->getParent().getDatasetId() == pointsDataset->getId();

                if (isValid)
                {
                    QVariant variant = _selectedClusterNamesVariant.getVariant();
                    QStringList clusterList = variant.toStringList();

                    QStringList speciesList;
                    auto speciesData = mv::data().getDataset<Clusters>(speciesDataset->getId());
                    auto species = speciesData->getClusters();
                    if (!species.empty())
                    {
                        for (auto& specie : species)
                        {
                            speciesList.push_back(specie.getName());
                        }
                    }
                    mv::Datasets filteredDatasets;
                    if (!speciesList.empty())
                    {
                        auto allDatasets = mv::data().getAllDatasets({ PointType });
                        for (const auto& dataset : allDatasets)
                        {
                            if (speciesList.contains(dataset->getGuiName()) && dataset->getDataType() == PointType)
                            {
                                filteredDatasets.push_back(dataset);
                            }
                        }
                    }

                    if (!filteredDatasets.empty())
                    {
                        for (const auto& dataset : filteredDatasets)
                        {
                            //qDebug() << dataset->getGuiName();
                            auto rawData = mv::data().getDataset < Points>(dataset.getDatasetId());
                            auto dimensionNames = rawData->getDimensionNames();
                            auto children = rawData->getChildren();
                            for (auto& child : children)
                            {
                                if (child->getGuiName() + "_mainData" == _hierarchyBottomClusterDataset.getCurrentDataset()->getGuiName())
                                {
                                    auto clustersData = mv::data().getDataset<Clusters>(child->getId());
                                    auto clusters = clustersData->getClusters();
                                    if (!clusters.empty())
                                    {
                                        std::vector<int> clusterIndicesSelected;
                                        std::vector<int> clusterIndicesFull(rawData->getNumPoints()); //fill it with rawData->getNumPoints()
                                        std::iota(clusterIndicesFull.begin(), clusterIndicesFull.end(), 0); // Fill the vector with increasing values
                                        for (auto& cluster : clusters)
                                        {
                                            if (clusterList.contains(cluster.getName()))
                                            {
                                                auto clusterIndices = cluster.getIndices();
                                                std::copy(clusterIndices.begin(), clusterIndices.end(), std::back_inserter(clusterIndicesSelected));
                                            }
                                        }
                                        if (!clusterIndicesSelected.empty())
                                        {
                                            for (int i = 0; i < dimensionNames.size(); i++)
                                            {
                                                std::vector<float> resultContainerShort(clusterIndicesSelected.size());
                                                std::vector<float> resultContainerFull(clusterIndicesFull.size());
                                                std::vector<int> dimensionIndex = { i };
                                                rawData->populateDataForDimensions(resultContainerShort, dimensionIndex, clusterIndicesSelected);
                                                rawData->populateDataForDimensions(resultContainerFull, dimensionIndex, clusterIndicesFull);
                                                float shortMean = calculateMean(resultContainerShort);
                                                float fullMean = calculateMean(resultContainerFull);
                                                _clusterNameToGeneNameToExpressionValue[dataset->getGuiName()][dimensionNames[i]] = shortMean / fullMean;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue());
            if (!geneListTable.isNull()) {
                _filteredGeneNamesVariant.setVariant(geneListTable);
            }
            else
            {
                qDebug()<<"No data found";
            }




        };
    connect(&_updateButtonForGeneFiltering, &TriggerAction::triggered, this, updateButtonForGeneFilteringUpdate);

    const auto speciesNamesDatasetUpdate = [this]() -> void
        {

        };
    connect(&_speciesNamesDataset, &DatasetPickerAction::currentIndexChanged, this, speciesNamesDatasetUpdate);

    const auto topNGenesFilterUpdate = [this]() -> void
        {

        };
    connect(&_topNGenesFilter, &IntegralAction::valueChanged, this, topNGenesFilterUpdate);


    _mainPointsDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyTopClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyMiddleClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyBottomClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _selectedClusterNamesVariant.setDefaultWidgetFlags(StringAction::WidgetFlag::LineEdit);
    _filteredGeneNamesVariant.setDefaultWidgetFlags(StringAction::WidgetFlag::LineEdit);
    _updateButtonForGeneFiltering.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);
    _speciesNamesDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _topNGenesFilter.setDefaultWidgetFlags(IntegralAction::WidgetFlag::SpinBox | IntegralAction::WidgetFlag::Slider);

}


SettingsAction::Widget::Widget(QWidget* parent, SettingsAction* SettingsAction) :
    WidgetActionWidget(parent, SettingsAction)
{ }

SettingsAction::OptionSelectionAction::Widget::Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction) :
    WidgetActionWidget(parent, optionSelectionAction)
{ }

inline SettingsAction::OptionSelectionAction::OptionSelectionAction(SettingsAction& SettingsAction) :
    GroupAction(nullptr, "CrossSpeciesComparisonGeneDetectPluginOptionSelectionAction"),
    _settingsAction(SettingsAction)
{
    setText("Options");
    setIcon(Application::getIconFont("FontAwesome").getIcon("wrench"));
    addAction(&_settingsAction.getMainPointsDataset());
    addAction(&_settingsAction.getHierarchyTopClusterDataset());
    addAction(&_settingsAction.getHierarchyMiddleClusterDataset());
    addAction(&_settingsAction.getHierarchyBottomClusterDataset());
    addAction(&_settingsAction.getSpeciesNamesDataset());
    addAction(&_settingsAction.getSelectedClusterNames());
    addAction(&_settingsAction.getFilteredGeneNames());
    addAction(&_settingsAction.getTopNGenesFilter());
    addAction(&_settingsAction.getUpdateButtonForGeneFiltering());
}


void SettingsAction::fromVariantMap(const QVariantMap& variantMap)
{
    WidgetAction::fromVariantMap(variantMap);

    _mainPointsDataset.fromParentVariantMap(variantMap);
    _hierarchyTopClusterDataset.fromParentVariantMap(variantMap);
    _hierarchyMiddleClusterDataset.fromParentVariantMap(variantMap);
    _hierarchyBottomClusterDataset.fromParentVariantMap(variantMap);
    _selectedClusterNamesVariant.fromParentVariantMap(variantMap);
    _filteredGeneNamesVariant.fromParentVariantMap(variantMap);
    //_updateButtonForGeneFiltering.fromParentVariantMap(variantMap);
    _topNGenesFilter.fromParentVariantMap(variantMap);
    _speciesNamesDataset.fromParentVariantMap(variantMap);
}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = WidgetAction::toVariantMap();

    _mainPointsDataset.insertIntoVariantMap(variantMap);
    _hierarchyTopClusterDataset.insertIntoVariantMap(variantMap);
    _hierarchyMiddleClusterDataset.insertIntoVariantMap(variantMap);
    _hierarchyBottomClusterDataset.insertIntoVariantMap(variantMap);
    _selectedClusterNamesVariant.insertIntoVariantMap(variantMap);
    _filteredGeneNamesVariant.insertIntoVariantMap(variantMap);
    //_updateButtonForGeneFiltering.insertIntoVariantMap(variantMap);
    _topNGenesFilter.insertIntoVariantMap(variantMap);
    _speciesNamesDataset.insertIntoVariantMap(variantMap);

    return variantMap;
}