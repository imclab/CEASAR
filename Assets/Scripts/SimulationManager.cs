﻿using System;
using System.Collections.Generic;
using System.Text;
using UnityEngine;

// This class is the manager Singleton, and contains specific references to application-level objects
public class SimulationManager
{
    // can't use constructor, guaranteed singleton
    protected SimulationManager() {
        LocalPlayer = new UserRecord();
        Debug.Log($"User has been created {LocalPlayer.Username}");
    }
    private static SimulationManager instance;

    public static SimulationManager GetInstance()
    {
        return instance ?? (instance = new SimulationManager());
    }

    // A decent random generator
    public System.Random rng = new System.Random();

    // List of scenes in project for current build target
    private readonly string[] _desktopScenes = new string[4] { "LoadSim", "Stars", "Horizon", "EarthInteraction" };
    private readonly string[] _ARScenes = new string[4] { "LoadSim", "Stars", "Horizon", "EarthInteraction" };
    private readonly string[] _VRScenes = new string[4] { "LoadSim", "Stars", "Horizon", "EarthInteraction" };

    public string[] Scenes {
        get {
#if UNITY_STANDALONE
            return this._desktopScenes;
#elif UNITY_IOS
            return this._ARScenes;
#elif UNITY_ANDROID
            return this._VRScenes;
#else
            // Catch-all default
            return this._desktopScenes;
#endif
        }
    }
    public GameObject NetworkControllerObject;
    public GameObject InteractionControllerObject;
    public ServerRecord server = null;
    private GameObject celestialSphereObject;

    public GameObject CelestialSphereObject
    {
        get { return celestialSphereObject; }
        set
        {
            celestialSphereObject = value;
            DataControllerComponent = celestialSphereObject.GetComponent<DataController>();
            MarkersControllerComponent = celestialSphereObject.GetComponentInChildren<MarkersController>();
            ConstellationsControllerComponent = celestialSphereObject.GetComponentInChildren<ConstellationsController>();
        }
    }

    public DataController DataControllerComponent { get; private set; }
    public MarkersController MarkersControllerComponent { get; private set; }
    public ConstellationsController ConstellationsControllerComponent { get; private set; }

    public bool IsReady = false;

    public string[] AnimalNames;
    public List<string> ColorNames = new List<string>();
    public List<Color> ColorValues = new List<Color>();

    public UserRecord LocalPlayer { get; set; }

    public Color LocalPlayerColor {
        get { return LocalPlayer.color; }
    }

    public string LocalUsername
    {
        get { return LocalPlayer.Username; }
        set {
            LocalPlayer = new UserRecord(value);
            Debug.Log($"User has been Changed: {LocalPlayer.Username}");
        }
    }
    public Color HorizonGroundColor = Color.green;

    // initial setup scale
    public readonly float InitialRadius = 100;
    // Scene Radius will be set in DataController, but we can keep a reference here for lookups elsewhere
    public float SceneRadius = 100;
    public float CurrentScaleFactor(float sceneRadius)
    {
        SceneRadius = sceneRadius;
        return sceneRadius / InitialRadius;
    }
    // Movement synchronization throttled for Heroku/Mongo
    public float MovementSendInterval = 1.0f;
    
    // Server Address!
    public string LocalNetworkServer = "ws://localhost:2567";
    public string DevNetworkServer = "ws://ceasar-serve-170822523-mftx7jq.herokuapp.com/";
    public string ProductionNetworkServer = "ws://ceasar-server-staging.herokuapp.com/";

    public StarComponent CurrentlySelectedStar;

    
    public float GetRelativeMagnitude(float starMagnitude)
    {
        //float min = DataControllerComponent.minMag;
        float max = DataControllerComponent.maxMag + Mathf.Abs(DataControllerComponent.minMag);
        return max - (starMagnitude + Mathf.Abs(DataControllerComponent.minMag));
    }
}