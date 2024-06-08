#pragma once
#include <actions/WidgetAction.h>
#include <actions/IntegralAction.h>
#include <actions/DecimalAction.h>
#include <actions/OptionAction.h>
#include <actions/OptionsAction.h>
#include <actions/ToggleAction.h>
#include "actions/DatasetPickerAction.h"
#include "PointData/PointData.h"
#include "ClusterData/ClusterData.h"
#include "event/EventListener.h"
#include "actions/Actions.h"
#include "Plugin.h"
#include "DataHierarchyItem.h"
#include "Set.h"
#include <AnalysisPlugin.h>
#include <memory>
#include <algorithm>    
#include <QDebug>
#include <QLabel>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include <QGridLayout>
#include <QFormLayout>
#include <QString>
#include <string>
#include <event/Event.h>
#include <QDebug>
#include <QLabel>
#include <string>
#include "actions/VariantAction.h"
#include "actions/GroupAction.h"
using namespace mv::gui;
class QMenu;
class CrossSpeciesComparisonClusterRankPlugin;

class FetchMetaData;
namespace mv
{
    class CoreInterface;
}

    class SettingsAction : public WidgetAction
    {
    public:
        class OptionSelectionAction : public GroupAction
        {
        protected:
            class Widget : public mv::gui::WidgetActionWidget {
            public:
                Widget(QWidget* parent, OptionSelectionAction* optionSelectionAction);
            };

            QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
                return new OptionSelectionAction::Widget(parent, this);
            };

        public:
            OptionSelectionAction(SettingsAction& SettingsAction);

        protected:
            SettingsAction& _settingsAction;
        };



    protected:

        class Widget : public mv::gui::WidgetActionWidget {
        public:
            Widget(QWidget* parent, SettingsAction* SettingsAction);
        };

        QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
            return new SettingsAction::Widget(parent, this);
        };

    public:
        SettingsAction(CrossSpeciesComparisonClusterRankPlugin& CrossSpeciesComparisonClusterRankPlugins);

public: // Action getters

    DatasetPickerAction& getMainPointsDataset() { return _mainPointsDataset; }
    DatasetPickerAction& getHierarchyTopClusterDataset() { return _hierarchyTopClusterDataset; }
    DatasetPickerAction& getHierarchyMiddleClusterDataset() { return _hierarchyMiddleClusterDataset; }
    DatasetPickerAction& getHierarchyBottomClusterDataset() { return _hierarchyBottomClusterDataset; }
    VariantAction& getSelectedClusterNames() { return _selectedClusterNamesVariant; }
    //VariantAction& getFilteredGeneNames() { return _filteredGeneNamesVariant; }
    //TriggerAction& getUpdateButtonForGeneFiltering() { return _updateButtonForGeneFiltering; }
    DatasetPickerAction& getSpeciesNamesDataset() { return _speciesNamesDataset; }
    DatasetPickerAction& getFilterTreeDataset() { return _filterTreeDataset; }
    //IntegralAction& getTopNGenesFilter() { return _topNGenesFilter; }
    OptionSelectionAction& getOptionSelectionAction() { return _optionSelectionAction; }
    DatasetPickerAction& getReferenceTreeDataset() { return _referenceTreeDataset; }
    StringAction& getGeneNamesConnection() { return _geneNamesConnection; }
    //DecimalAction& getTreeSimilarity() { return _treeSimilarity; }
    double* condensedDistanceMatrix(std::vector<float>& items);
    std::string mergeToNewick(int* merge, int numOfLeaves);
    QString createJsonTreeFromNewick(QString tree, std::vector<QString> leafNames);
    private:
        QVariant createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map, const QString& treeDatasetId, const float& treeSimilarityScore, const std::map<QString, std::vector<QString>>& geneCounter, const int& n);
        QVariant findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n, QString datasetId, float treeSimilarityScore);

public: // Serialization

    /**
     * Load widget action from variant map
     * @param Variant map representation of the widget action
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save widget action to variant map
     * @return Variant map representation of the widget action
     */
    QVariantMap toVariantMap() const override;

protected:
    CrossSpeciesComparisonClusterRankPlugin& _crossSpeciesComparisonClusterRankPlugin;
    DatasetPickerAction    _mainPointsDataset;
    DatasetPickerAction    _hierarchyTopClusterDataset;
    DatasetPickerAction    _hierarchyMiddleClusterDataset;
    DatasetPickerAction    _hierarchyBottomClusterDataset;
    VariantAction           _selectedClusterNamesVariant;
    //VariantAction           _filteredGeneNamesVariant;
    //TriggerAction          _updateButtonForGeneFiltering;
    DatasetPickerAction    _speciesNamesDataset;
    DatasetPickerAction    _filterTreeDataset;
    //std::map<QString, std::map<QString, float>> _clusterNameToGeneNameToExpressionValue;
    //IntegralAction          _topNGenesFilter;
    OptionSelectionAction         _optionSelectionAction;
    DatasetPickerAction    _referenceTreeDataset;
    StringAction           _geneNamesConnection;
    //DecimalAction          _treeSimilarity;

};