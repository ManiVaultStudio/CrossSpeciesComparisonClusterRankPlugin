#include "SettingsAction.h"
#include <vector>
#include <numeric>
#include <execution>
#include <iostream>
#include <algorithm>
#include <vector>


void printMap(const std::map<QString, std::map<QString, float>>& map) {
    for (const auto& outerPair : map) {
        std::cout << "Cluster Name: " << outerPair.first.toStdString() << std::endl;
        for (const auto& innerPair : outerPair.second) {
            std::cout << "    Gene Name: " << innerPair.first.toStdString() << ", Expression Value: " << innerPair.second << std::endl;
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
SettingsAction::SettingsAction(QObject* parent) :
    GroupAction(parent, "SettingsAction", true),
    _mainPointsDataset(this, "Main Points Dataset"),
    _hierarchyTopClusterDataset(this, "Hierarchy Top Cluster Dataset"),
    _hierarchyMiddleClusterDataset(this, "Hierarchy Middle Cluster Dataset"),
    _hierarchyBottomClusterDataset(this, "Hierarchy Bottom Cluster Dataset"),
    _selectedClusterNamesVariant(this, "Selected Cluster Names Variant"),
    _filteredGeneNamesVariant(this, "Filtered Gene Names Variant"),
    _updateButtonForGeneFiltering(this, "Filter Genes Trigger"),
    _speciesNamesDataset(this, "Species Names Dataset")
{
    setText("Cross-Species Comparison Cluster Rank Settings");
    _mainPointsDataset.setToolTip("Main Points Dataset");
    _hierarchyTopClusterDataset.setToolTip("Hierarchy Top Cluster Dataset");
    _hierarchyMiddleClusterDataset.setToolTip("Hierarchy Middle Cluster Dataset");
    _hierarchyBottomClusterDataset.setToolTip("Hierarchy Bottom Cluster Dataset");
    _selectedClusterNamesVariant.setToolTip("Selected Cluster Names Variant");
    _filteredGeneNamesVariant.setToolTip("Filtered Gene Names Variant");
    _updateButtonForGeneFiltering.setToolTip("Filter Genes Trigger");
    _speciesNamesDataset.setToolTip("Species Names Dataset");



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
            auto speciesDataset= _speciesNamesDataset.getCurrentDataset();
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
                    if (species.size()>0)
                    {
                        for (auto specie : species)
                        {
                            speciesList.push_back(specie.getName());
                        }
                    }
                    mv::Datasets filteredDatasets;
                    if (speciesList.size() > 0)
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


                    if (filteredDatasets.size()>0)
                    {
                    for (const auto& dataset : filteredDatasets)
                    {
                        qDebug() << dataset->getGuiName();
                        auto rawData = mv::data().getDataset < Points>(dataset.getDatasetId());
                        auto dimensionNames= rawData->getDimensionNames();
                        auto children = rawData->getChildren();
                        for (auto child : children)
                        {
                            if (child->getGuiName() + "_mainData" == _hierarchyBottomClusterDataset.getCurrentDataset()->getGuiName())
                            {
                                auto clustersData = mv::data().getDataset<Clusters>(child->getId());
                                auto clusters = clustersData->getClusters();
                                if (clusters.size() > 0)
                                {
                                    std::vector<int> clusterIndicesSelected;
                                    
                                    std::vector<int> clusterIndicesFull;//fill it with rawData->getNumPoints()
                                    for (int i = 0; i < rawData->getNumPoints(); i++)
                                    {
                                        clusterIndicesFull.push_back(i);
                                    }
                                    for (auto cluster : clusters)
                                    {
                                        if (clusterList.contains(cluster.getName()) )
                                        {
                                            auto clusterIndices = cluster.getIndices();
                                            std::copy(clusterIndices.begin(), clusterIndices.end(), std::back_inserter(clusterIndicesSelected));          
                                        }

                                    }
                                    if (clusterIndicesSelected.size()>0)
                                    {
                                        for (int i = 0; i < dimensionNames.size(); i++)
                                        {
                                            std::vector<float> resultContainerShort(clusterIndicesSelected.size() * 1);
                                            std::vector<float> resultContainerFull(clusterIndicesFull.size() * 1);
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

            printMap(_clusterNameToGeneNameToExpressionValue);


        };
    connect(&_updateButtonForGeneFiltering, &TriggerAction::triggered, this, updateButtonForGeneFilteringUpdate);

    const auto speciesNamesDatasetUpdate = [this]() -> void
        {

        };
    connect(&_speciesNamesDataset, &DatasetPickerAction::currentIndexChanged, this, speciesNamesDatasetUpdate);

    _mainPointsDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyTopClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyMiddleClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _hierarchyBottomClusterDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);
    _selectedClusterNamesVariant.setDefaultWidgetFlags(StringAction::WidgetFlag::LineEdit);
    _filteredGeneNamesVariant.setDefaultWidgetFlags(StringAction::WidgetFlag::LineEdit);
    _updateButtonForGeneFiltering.setDefaultWidgetFlags(TriggerAction::WidgetFlag::IconText);
    _speciesNamesDataset.setDefaultWidgetFlags(DatasetPickerAction::WidgetFlag::ComboBox);

    setSerializationName("CSEC:Cross-Species Comparison Cluster Rank Settings");

    _mainPointsDataset.setSerializationName("CSCCR:MainPointsDataset");
    _hierarchyTopClusterDataset.setSerializationName("CSCCR:HierarchyTopClusterDataset");
    _hierarchyMiddleClusterDataset.setSerializationName("CSCCR:HierarchyMiddleClusterDataset");
    _hierarchyBottomClusterDataset.setSerializationName("CSCCR:HierarchyBottomClusterDataset");
    _selectedClusterNamesVariant.setSerializationName("CSCCR:SelectedClusterNamesVariant");
    _filteredGeneNamesVariant.setSerializationName("CSCCR:FilteredGeneNamesVariant");
    //_updateButtonForGeneFiltering.setSerializationName("CSCCR:UpdateButtonForGeneFiltering");
    _speciesNamesDataset.setSerializationName("CSCCR:SpeciesNamesDataset");

    addAction(&_mainPointsDataset);
    addAction(&_hierarchyTopClusterDataset);
    addAction(&_hierarchyMiddleClusterDataset);
    addAction(&_hierarchyBottomClusterDataset);
    addAction(&_speciesNamesDataset);
    addAction(&_selectedClusterNamesVariant);
    addAction(&_filteredGeneNamesVariant);
    addAction(&_updateButtonForGeneFiltering);




}

void SettingsAction::fromVariantMap(const QVariantMap& variantMap)
{
    GroupAction::fromVariantMap(variantMap);


    _mainPointsDataset.fromParentVariantMap(variantMap);
    _hierarchyTopClusterDataset.fromParentVariantMap(variantMap);
    _hierarchyMiddleClusterDataset.fromParentVariantMap(variantMap);
    _hierarchyBottomClusterDataset.fromParentVariantMap(variantMap);
    _selectedClusterNamesVariant.fromParentVariantMap(variantMap);
    _filteredGeneNamesVariant.fromParentVariantMap(variantMap);
    //_updateButtonForGeneFiltering.fromParentVariantMap(variantMap);
    _speciesNamesDataset.fromParentVariantMap(variantMap);

}

QVariantMap SettingsAction::toVariantMap() const
{
    QVariantMap variantMap = GroupAction::toVariantMap();


    _mainPointsDataset.insertIntoVariantMap(variantMap);
    _hierarchyTopClusterDataset.insertIntoVariantMap(variantMap);
    _hierarchyMiddleClusterDataset.insertIntoVariantMap(variantMap);
    _hierarchyBottomClusterDataset.insertIntoVariantMap(variantMap);
    _selectedClusterNamesVariant.insertIntoVariantMap(variantMap);
    _filteredGeneNamesVariant.insertIntoVariantMap(variantMap);
    //_updateButtonForGeneFiltering.insertIntoVariantMap(variantMap);
    _speciesNamesDataset.insertIntoVariantMap(variantMap);


    return variantMap;
}