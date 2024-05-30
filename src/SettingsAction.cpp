#include "SettingsAction.h"
#include "CrossSpeciesComparisonClusterRankPlugin.h"
#include <vector>
#include <numeric>
#include <execution>
#include <iostream>
#include <algorithm>
#include <vector>
#include <CrossSpeciesComparisonTreeData.h>
#include "lib/Distance/annoylib.h"
#include "lib/Distance/kissrandom.h"
#include <QtConcurrent>
#include "lib/JSONnlohmann/json.hpp"
#include "lib/Clustering/fastcluster.h"
#include <stack>
#include <sstream>
float calculateVariance(const std::vector<float>& numbers) {
    float sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
    float mean = sum / numbers.size();
    float variance = 0.0;

    for (const auto& num : numbers) {
        variance += std::pow(num - mean, 2);
    }

    return variance / numbers.size();
}

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
struct Statistics {
    float mean; 
    float variance;
    float stdDeviation;
};

Statistics calculateStatistics(const std::vector<float>& numbers) {
    if (numbers.empty()) {
        return { 0.0, 0.0, 0.0 };
    }

    float sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
    float mean = sum / numbers.size();
    float sq_sum = std::inner_product(numbers.begin(), numbers.end(), numbers.begin(), 0.0);
    float variance = sq_sum / numbers.size() - mean * mean;
    float stdDeviation = std::sqrt(variance);

    return { mean, variance, stdDeviation };
}

std::string mergeToNewick(int* merge, int numOfLeaves) {
    std::vector<std::string> labels(numOfLeaves);
    for (int i = 0; i < numOfLeaves; ++i) {
        labels[i] = std::to_string(i + 1);
    }

    std::stack<std::string> stack;

    for (int i = 0; i < 2 * (numOfLeaves - 1); i += 2) {
        int left = merge[i];
        int right = merge[i + 1];

        std::string leftStr;
        if (left < 0) {
            leftStr = labels[-left - 1];
        }
        else {
            leftStr = stack.top();
            stack.pop();
        }

        std::string rightStr;
        if (right < 0) {
            rightStr = labels[-right - 1];
        }
        else {
            rightStr = stack.top();
            stack.pop();
        }

        std::string merged = "(" + leftStr + "," + rightStr + ")";
        stack.push(merged);
    }

    return stack.top() + ";";
}
double* condensedDistanceMatrix(std::vector<float>& items) {
    int n = items.size();
    double* distmat = new double[(n * (n - 1)) / 2];
    int k = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            distmat[k] = std::abs(items[i] - items[j]);
            ++k;
        }
    }
    return distmat;
   
    //    int n = items.size();
    //double* distmat = new double[(n * (n - 1)) / 2];
    //int k = 0;

    //// Build Annoy index
    //Annoy::AnnoyIndex<int32_t, float, Annoy::Euclidean, Annoy::Kiss32Random, Annoy::AnnoyIndexSingleThreadedBuildPolicy> index(1);
    ////Annoy::AnnoyIndex<int32_t, float, Annoy::Manhattan, Annoy::Kiss32Random, Annoy::AnnoyIndexSingleThreadedBuildPolicy> index(1);
    ////Annoy::AnnoyIndex<int32_t, float, Annoy::Angular, Annoy::Kiss32Random, Annoy::AnnoyIndexSingleThreadedBuildPolicy> index(1);
    //for (int i = 0; i < n; ++i) {
    //    index.add_item(i, &items[i]);
    //}
    //index.build(10);  // 10 is the number of trees for the index. More trees gives higher precision.

    //// Calculate distances
    //for (int i = 0; i < n; ++i) {
    //    for (int j = i + 1; j < n; ++j) {
    //        distmat[k] = index.get_distance(i, j);
    //        ++k;
    //    }
    //}

    //return distmat;
 
}

QVariant createModelFromData(const QStringList& returnGeneList, const std::map<QString, std::map<QString, float>>& map, std::vector<QString> leafnames) {
    
    if (returnGeneList.isEmpty() || map.empty()) {
        return QVariant();
    }

    QStandardItemModel* model = new QStandardItemModel();
    QStringList initColumnNames = { "ID",  "Variance","Standard Deviation","Grand Mean"};
    model->setHorizontalHeaderLabels(initColumnNames);
    int numOfSpecies = map.size();
    std::map<QString, std::map<QString, float>>::const_iterator it = map.begin();
    for (int i = 0 + initColumnNames.size(); i < numOfSpecies + initColumnNames.size(); i++, it++) {
        QString headerTitle =  it->first;
        headerTitle.replace("_", " ");
        headerTitle = QString("Mean_") + headerTitle;
        model->setHorizontalHeaderItem(i, new QStandardItem(headerTitle));
    }
    std::map<QString, QString> newickTrees;
    for (auto gene : returnGeneList)
    {
        QList<QStandardItem*> row;
        std::vector<float> numbers;


        for (const auto& outerPair : map) {
            QString outerKey = outerPair.first;
            const std::map<QString, float>& innerMap = outerPair.second;
            try {
                numbers.push_back(innerMap.at(gene));
            }
            catch (const std::out_of_range& e) {
                numbers.push_back(0);
            }

            

        }

        const std::unordered_map<std::string, int> clusteringTypeMap = {
    {"Complete", HCLUST_METHOD_COMPLETE},
    {"Average", HCLUST_METHOD_AVERAGE},
    {"Median", HCLUST_METHOD_MEDIAN}
        };
        std::string clusteringTypecurrentText = "Single";  //"Single","Complete", "Average","Median"
        int opt_method = clusteringTypeMap.count(clusteringTypecurrentText) ? clusteringTypeMap.at(clusteringTypecurrentText) : HCLUST_METHOD_SINGLE;

        auto numOfLeaves = numOfSpecies;
        double* distmat = new double[(numOfLeaves * (numOfLeaves - 1)) / 2];

        
        distmat = condensedDistanceMatrix(numbers);

        int* merge = new int[2 * (numOfLeaves - 1)];
        double* height = new double[numOfLeaves - 1];
        hclust_fast(numOfLeaves, distmat, opt_method, merge, height);
        std::string newick = mergeToNewick(merge, numOfLeaves);
        int totalChars = newick.length();
        //std::cout << "\nOriginal Newick format: " << newick << std::endl;
        //add gene and newick to newickTrees
        newickTrees.insert(std::make_pair(gene, QString::fromStdString(newick)));

        /*
        int i = 0;
        newick += ';';
        std::string jsonString = "";
        std::stringstream jsonStream;
        while (i < newick.size()) {
            if (newick[i] == '(') {
                jsonStream << "{\n\"children\": [";
                i++;
            }
            else if (newick[i] == ',') {
                jsonStream << ",";
                i++;
            }
            else if (newick[i] == ')') {
                jsonStream << "],\n\"id\": 1,\n\"score\": 1,\n\"width\": 1\n}";
                i++;
            }
            else if (newick[i] == ';') {
                break;
            }
            else {
                if (isdigit(newick[i])) {
                    int skip = 1;
                    std::string num = "";
                    for (int j = i; j < newick.size(); j++) {
                        if (isdigit(newick[j])) {
                            continue;
                        }
                        else {
                            num = newick.substr(i, j - i);

                            skip = j - i;
                            break;
                        }
                    }


                    std::string species = leafnames[(std::stoi(num) - 1)].toStdString();
                    jsonStream << "{\n\"color\": \"#000000\",\n\"hastrait\": true,\n\"iscollapsed\": false,\n\"name\": \"" << species << "\"\n}";
                    i += skip;
                }
            }
        }

        jsonString = jsonStream.str();

        nlohmann::json json = nlohmann::json::parse(jsonString);
        std::string jsonStr = json.dump(4);

        QJsonObject valueStringReference = QJsonDocument::fromJson(QString::fromStdString(jsonStr).toUtf8()).object();

        QString completedString = QJsonDocument(valueStringReference).toJson(QJsonDocument::Compact);

        */
delete[] distmat;
delete[] merge;
delete[] height;

Statistics stats = calculateStatistics(numbers);

row.push_back(new QStandardItem(gene));
row.push_back(new QStandardItem(QString::number(stats.variance)));
row.push_back(new QStandardItem(QString::number(stats.stdDeviation)));
row.push_back(new QStandardItem(QString::number(stats.mean)));
for (auto numb : numbers)
{
    row.push_back(new QStandardItem(QString::number(numb)));
}



model->appendRow(row);
    }




    //check which newick trees are the same and group them together
    std::map<QString, std::vector<QString>> clusteringMap;




    for (auto it = newickTrees.begin(); it != newickTrees.end(); ++it) {
        QString currentNewick = it->second;
        QString currentGene = it->first;
        if (clusteringMap.empty()) {
            clusteringMap[currentNewick] = { currentGene };
        }
        else {
            bool found = false;
            for (auto& cluster : clusteringMap) {
                QString clusterNewick = cluster.first;
                if (currentNewick == clusterNewick) {
                    cluster.second.push_back(currentGene);
                    found = true;
                    break;
                }
            }
            if (!found) {
                clusteringMap[currentNewick] = { currentGene };
            }
        }
    }
    //add the genes ABLIM3, ADCK2, ADCY10 in the vector std::vector<QString> to the first key of clusteringMap for testing
    //clusteringMap[clusteringMap.begin()->first] = { "ABLIM3", "ADCK2", "ADCY10" };
    //print clusteringMap    '
    for (auto& cluster : clusteringMap) {
        QString newick = cluster.first;
        std::vector<QString> genes = cluster.second;

        if (genes.size() > 1) {
            std::cout << "Newick: " << newick.toStdString() << std::endl;
            std::cout << "Genes: ";
            for (auto& gene : genes) {
                std::cout << gene.toStdString() << ", ";
            }
            std::cout << std::endl;
        }

    }





    //color the rows of the genes that have the same newick tree using color codes  #8dd3c7#ffffb3#bebada#fb8072#80b1d3#fdb462#b3de69#fccde5#d9d9d9#bc80bd#ccebc5#ffed6f and    #a6cee3#1f78b4#b2df8a#33a02c#fb9a99#e31a1c#fdbf6f#ff7f00#cab2d6#6a3d9a#ffff99#b15928
    QStringList colorCodes = { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928" };
    int colorIndex = 0;
    for (auto& cluster : clusteringMap) {
        QString newick = cluster.first;
        std::vector<QString> genes = cluster.second;

        if (genes.size() > 1) {
            for (auto& gene : genes) {
                for (int i = 0; i < model->rowCount(); i++) {
                    if (model->item(i, 0)->text() == gene) {
                        model->item(i, 0)->setBackground(QBrush(QColor(colorCodes[colorIndex])));
                    }
                }
            }
            colorIndex++;
        }
    }



    return QVariant::fromValue(model);


    //qRegisterMetaType<QStandardItemModel>("QStandardItemModel");
    //QVariant variant = QVariant::fromValue(model);
    //return variant;
}

QVariant  findTopNGenesPerCluster(const std::map<QString, std::map<QString, float>>& map, int n, std::vector<QString> leafnames) {
    
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
    return createModelFromData(returnGeneList, map, leafnames);
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
            std::vector<QString> leafnames;
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
                            leafnames.push_back(dataset->getGuiName());
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
                                                float meanValue = 0.0;
                                                if (fullMean != 0.0)
                                                {
                                                    meanValue = shortMean / fullMean;
                                                }
                                                _clusterNameToGeneNameToExpressionValue[dataset->getGuiName()][dimensionNames[i]] = meanValue;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            QVariant geneListTable = findTopNGenesPerCluster(_clusterNameToGeneNameToExpressionValue, _topNGenesFilter.getValue(), leafnames);
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