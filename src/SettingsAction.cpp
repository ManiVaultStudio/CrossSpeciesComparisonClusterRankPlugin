#include "SettingsAction.h"



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
            
            if (clusterDataset.isValid() && pointsDataset.isValid() && speciesDataset.isValid())
            {
                isValid = clusterDataset->getParent().getDatasetId() == pointsDataset->getId() && speciesDataset->getParent().getDatasetId() == pointsDataset->getId();

                if (isValid)
                {
                    QVariant variant = _selectedClusterNamesVariant.getVariant();
                    QStringList stringList = variant.toStringList();
                    qDebug() << stringList.join(", ");
                }
            }



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