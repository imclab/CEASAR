﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;
using System.Text;

public class StarInfoPanel : MonoBehaviour
{
    private SimulationManager manager;
    public TextMeshProUGUI starInfoText;

    public void Setup(MainUIController mainUIController)
    {
        manager = SimulationManager.GetInstance();
        if (manager.CurrentlySelectedStar == null)
        {
            starInfoText.text = "";
            mainUIController.HidePanel("StarInfoPanel");
        }
        else if (string.IsNullOrEmpty(starInfoText.text))
        {
            UpdateStarInfoPanel();
            Debug.Log("Highlighting constellation " + manager.CurrentlySelectedStar.starData.ConstellationFullName);
            ConstellationsController constellationsController = FindObjectOfType<ConstellationsController>();
            if (constellationsController)
            {
                constellationsController.HighlightSingleConstellation(manager.CurrentlySelectedStar.starData.ConstellationFullName);
            }
        }
    }

    public void UpdateStarInfoPanel()
    {
        manager = SimulationManager.GetInstance();
        if (manager.CurrentlySelectedStar != null)
        {
            Star starData = manager.CurrentlySelectedStar.starData;
            StringBuilder description = new StringBuilder();
            description.Append("Name: ").Append(starData.ProperName.Length > 0 ? starData.ProperName : "N/A");
            description.Append("   Constellation: ").AppendLine(starData.ConstellationFullName.Length > 0 ? starData.ConstellationFullName : "N/A");
            description.Append("Hipparcos #: ").AppendLine(starData.Hipparcos.ToString());
            description.Append("Bayer Des: ").AppendLine(starData.BayerDesignation.Length > 0 ? starData.BayerDesignation : "N/A");
            description.Append("Flamsteed Des: ").AppendLine(starData.FlamsteedDesignation.Length > 0 ? starData.FlamsteedDesignation : "N/A");
            description.Append("m: ").AppendLine(starData.Mag.ToString());
            starInfoText.text = description.ToString();
        }
    }
    public void ClearStarSelection()
    {
        manager.CurrentlySelectedStar = null;
        starInfoText.text = "";
        ConstellationsController constellationsController = FindObjectOfType<ConstellationsController>();
        if (constellationsController)
        {
            // TODO: We network selecting stars, we need to network de-selecting stars also
            constellationsController.HighlightAllConstellations(true);
        }
    }
}
