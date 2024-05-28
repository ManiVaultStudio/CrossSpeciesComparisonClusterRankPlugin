#pragma once

#include <actions/DatasetPickerAction.h>
#include "PointData/DimensionsPickerAction.h"
#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include "actions/GroupAction.h"
#include "actions/TriggerAction.h"
#include "actions/ToggleAction.h"
#include "actions/OptionsAction.h"
#include "actions/OptionAction.h"
#include "actions/StringAction.h"
#include "actions/VariantAction.h"

/** All GUI related classes are in the ManiVault Graphical User Interface namespace */
using namespace mv::gui;

/**
 * CrossSpeciesComparisonSettingsAction class
 *
 * Class that houses settings for the CrossSpeciesComparison analysis plugin
 *
 * This settings class is derived from the group action class. A group action
 * is a special type of action; when injected into a dataset its user interface
 * becomes available in the data properties widget. The group action list the child
 * actions in a form-like fashion. The order in which they appear corresponds with
 * the order of declaration.
 *
 * Note: we strongly encourage you to use ManiVault core actions to build the user
 * interface. Actions separate the data and business logic from the user interface.
 * We have standard actions for editing of strings, decimals, integrals, options,
 * color and color maps. With these components, there is no need to write to create
 * the user interface yourself. The actions will take care of this. For more
 * information regarding actions, please visit actions/Actions.h
 */
class SettingsAction : public GroupAction
{
public:

    /**
     * Constructor
     * @param parent Pointer to parent object
     */
    SettingsAction(QObject* parent = nullptr);

public: // Action getters


    DatasetPickerAction& getMainPointsDataset() { return _mainPointsDataset; }
    DatasetPickerAction& getHierarchyTopClusterDataset() { return _hierarchyTopClusterDataset; }
    DatasetPickerAction& getHierarchyMiddleClusterDataset() { return _hierarchyMiddleClusterDataset; }
    DatasetPickerAction& getHierarchyBottomClusterDataset() { return _hierarchyBottomClusterDataset; }
    VariantAction& getSelectedClusterNames() { return _selectedClusterNamesVariant; }
    VariantAction& getFilteredGeneNames() { return _filteredGeneNamesVariant; }
    TriggerAction& getUpdateButtonForGeneFiltering() { return _updateButtonForGeneFiltering; }
    DatasetPickerAction& getSpeciesNamesDataset() { return _speciesNamesDataset; }


    
    

public: // Serialization

    /**
     * Load plugin from variant map
     * @param Variant map representation of the plugin
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save plugin to variant map
     * @return Variant map representation of the plugin
     */
    QVariantMap toVariantMap() const override;

public:
    DatasetPickerAction    _mainPointsDataset;
    DatasetPickerAction    _hierarchyTopClusterDataset;
    DatasetPickerAction    _hierarchyMiddleClusterDataset;
    DatasetPickerAction    _hierarchyBottomClusterDataset;
    VariantAction           _selectedClusterNamesVariant;
    VariantAction           _filteredGeneNamesVariant;
    TriggerAction          _updateButtonForGeneFiltering;
    DatasetPickerAction    _speciesNamesDataset;
    std::map<QString, std::map<QString, float>> _clusterNameToGeneNameToExpressionValue;

};
