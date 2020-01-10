﻿using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using System.Collections.Generic;
using GameDevWare.Serialization;
using UnityEngine.SceneManagement;




public enum NetworkConnection { Local, Dev, Remote, None }

[RequireComponent(typeof(ColyseusClient))]
public class NetworkController : MonoBehaviour
{
    public const string PLAYER_PREFS_NAME_KEY = "CAESAR_USERNAME";
    // UI Buttons are attached through Unity Inspector
    public Button connectButton;
    public Button randomizeUsernameButton;
    public TMPro.TMP_Text connectButtonText;
    public InputField m_EndpointField;
    public Text connectionStatusText;
    public Text usernameText;
    private string _connectionStatusMessage;
    private NetworkConnection _selectedNetwork = NetworkConnection.None;

    public GameObject avatar;
    private Player localPlayer;
    private GameObject localPlayerAvatar;
    public GameObject interactionIndicator;

    // todo: set up 4-digit pins for a collab room
    public string roomName = "ceasar";

    protected ColyseusClient colyseusClient;

    protected IndexedDictionary<string, GameObject> remotePlayers = new IndexedDictionary<string, GameObject>();

    public TMPro.TMP_Text debugMessages;

    private string localUsername = "";

    public bool autoConnect = false;
    public List<string> scenesWithAvatars;

    public NetworkConnection networkConnection = NetworkConnection.Remote;

    public bool IsConnected
    {
        get { return colyseusClient != null && colyseusClient.IsConnected; }
    }

    public string ServerStatusMessage
    {
        set
        {
            _connectionStatusMessage = value;
            if (connectionStatusText) connectionStatusText.text = value;
            if (connectButtonText && IsConnected) connectButtonText.text = "Disconnect";
            if (connectButtonText && !IsConnected) connectButtonText.text = "Connect";
        }
    }
    public void NetworkPanelToggled(bool active)
    {
        if (active)
        {
            // refresh connection status
            ServerStatusMessage = _connectionStatusMessage;
            debugMessages.text = colyseusClient.GetClientList();
        }
    }

    // Called via event handlers set up in the scene in the NetworkController prefab
    public void SetNetworkAddress(string destination)
    {
        switch (destination)
        {
            case "local":
                setNetworkAddress(NetworkConnection.Local);
                break;
            case "dev":
                setNetworkAddress(NetworkConnection.Dev);
                break;
            default:
                setNetworkAddress(NetworkConnection.Remote);
                break;
        }
    }
    private void setNetworkAddress(NetworkConnection destination)
    {
        SimulationManager manager = SimulationManager.GetInstance();
        switch (destination)
        {
            case NetworkConnection.Local:
                m_EndpointField.text = manager.LocalNetworkServer;
                break;
            case NetworkConnection.Dev:
                m_EndpointField.text = manager.DevNetworkServer;
                break;
            case NetworkConnection.Remote:
                m_EndpointField.text = manager.ProductionNetworkServer;
                break;
            default:
                m_EndpointField.text = manager.ProductionNetworkServer;
                break;
        }

        _selectedNetwork = destination;
    }

    // Need this so the network UI persists across scenes
    private void Awake()
    {
        DontDestroyOnLoad(this.gameObject);
        EnsureUsername();
    }

    void EnsureUsername()
    {
        if(!string.IsNullOrEmpty(localUsername))
        {
            return;
        }
        string foundName = PlayerPrefs.GetString(PLAYER_PREFS_NAME_KEY);
        if (string.IsNullOrEmpty(foundName))
        {
            RandomizeUsername();
        }
        else
        {
            SetUsername(foundName);
        }
    }

    void Start()
    {
        colyseusClient = GetComponent<ColyseusClient>();
        /* Demo UI */
        connectButton.onClick.AddListener(ConnectToServer);

        if (autoConnect)
        {
            ConnectToServer();
        }
        SceneManager.sceneLoaded += OnSceneLoaded;
    }

    void OnSceneLoaded(Scene scene, LoadSceneMode mode)
    {
        Debug.Log("OnSceneLoaded: " + scene.name);
        EnsureUsername();
        if (IsConnected && scenesWithAvatars.Contains(scene.name))
        {
            OnPlayerAdd(localPlayer);
            foreach (string p in remotePlayers.Keys)
            {
                Player remotePlayer = colyseusClient.GetPlayerById(p);
                OnPlayerAdd(remotePlayer);
            }
        }
        CCLogger.Log(CCLogger.EVENT_SCENE, "OnSceneLoaded: " + scene.name);
    }

    void SetUsername(string uName)
    {
        SimulationManager manager = SimulationManager.GetInstance();
        localUsername = uName;
        manager.LocalPlayerColor = manager.GetColorForUsername(localUsername);
        localPlayerAvatar = GameObject.FindWithTag("LocalPlayerAvatar");
        UpdateLocalPlayerAvatar(uName);
        Debug.Log(localUsername);
        PlayerPrefs.SetString(PLAYER_PREFS_NAME_KEY, localUsername);
        if (usernameText)
        {
            usernameText.text = localUsername;
            usernameText.color = manager.LocalPlayerColor;
        }
        CCLogger.Log(CCLogger.EVENT_USERNAME, "Username set " + localUsername);
    }

    public void RandomizeUsername()
    {
        SimulationManager manager = SimulationManager.GetInstance();
        SetUsername(manager.GenerateUsername());
    }

    void ConnectToServer()
    {
        /*
         * Get Colyseus endpoint from InputField
         *  for localhost use ws://localhost:2567/        
         */
        if (!IsConnected)
        {
            SimulationManager manager = SimulationManager.GetInstance();
            EnsureUsername();
            string _localEndpoint = manager.LocalNetworkServer;
            string _remoteEndpoint = manager.ProductionNetworkServer;
            string endpoint = string.Empty;

            if (string.IsNullOrEmpty(m_EndpointField.text))
            {
                // no user interaction with the network address, work on some defaults
#if UNITY_EDITOR
                endpoint = _localEndpoint;
                if (_selectedNetwork == NetworkConnection.None) _selectedNetwork = NetworkConnection.Local;
#else
                endpoint = _remoteEndpoint;
                if (_selectedNetwork == NetworkConnection.None) _selectedNetwork = NetworkConnection.Remote;
#endif
            }
            else
            {
                endpoint = m_EndpointField.text;
            }
            colyseusClient.ConnectToServer(endpoint, localUsername);
            manager.NetworkStatus = _selectedNetwork;
            CCLogger.Log(CCLogger.EVENT_CONNECT, "connected");
            if (randomizeUsernameButton != null)  randomizeUsernameButton.enabled = false;
        }
        else if (IsConnected)
        {
            Disconnect();
        }
        else
        {
            Debug.Log("Already disconnected");
        }
    }

    void Disconnect()
    {
        colyseusClient.Disconnect();
        randomizeUsernameButton.enabled = true;
        // Destroy player game objects
        foreach (KeyValuePair<string, GameObject> entry in remotePlayers)
        {
            Destroy(entry.Value);
        }
        remotePlayers.Clear();
        CCLogger.Log(CCLogger.EVENT_DISCONNECT, "disconnected");
        debugMessages.text = "";
    }

    void UpdateLocalPlayerAvatar(string username)
    {
        if (localPlayerAvatar == null)
        {
            localPlayerAvatar = GameObject.FindWithTag("LocalPlayerAvatar");
        }
        Color playerColor = SimulationManager.GetInstance().LocalPlayerColor;
        if (localPlayerAvatar)
        {
            localPlayerAvatar.GetComponent<Renderer>().material.color = playerColor;
            localPlayerAvatar.name = "localPlayer_" + username;
        }
    }

    void updatePlayerList()
    {
        debugMessages.text = colyseusClient.GetClientList();
        debugMessages.enabled = true;
        string listOfPlayersForDebug = "";
        foreach (var p in remotePlayers.Keys)
        {
            listOfPlayersForDebug = listOfPlayersForDebug + p + " \n";
        }
        Debug.Log(listOfPlayersForDebug);
        debugMessages.text = listOfPlayersForDebug;
    }

    public void OnPlayerAdd(Player player)
    {
        bool isLocal = localUsername == player.username;
        Debug.Log("Player add! playerName: " + player.username + " playerId: " + player.id);
        Debug.Log("is local: " + isLocal);
        Debug.Log("localPlayer: " + localPlayer?.username);
        Vector3 pos = Utils.NetworkPosToPosition(player.playerPosition.position);
        Quaternion rot = Utils.NetworkRotToRotation(player.playerPosition.rotation);
        if (isLocal)
        {
            localPlayer = player;

            UpdateLocalPlayerAvatar(player.username);
        }
        else
        {
            GameObject remotePlayerAvatar = Instantiate(avatar, pos, rot);
            Color playerColor = SimulationManager.GetInstance().GetColorForUsername(player.username);
            remotePlayerAvatar.GetComponent<Renderer>().material.color = playerColor;
            remotePlayerAvatar.name = "remotePlayer_" + player.username;
            remotePlayerAvatar.AddComponent<RemotePlayerMovement>();

            // add "player" to map of players
            if (!remotePlayers.ContainsKey(player.username))
            {
                remotePlayers.Add(player.username, remotePlayerAvatar);
            }
            else
            {
                remotePlayers[player.username] = remotePlayerAvatar;
            }
        }
        updatePlayerList();
    }

    public void OnPlayerRemove(Player player)
    {
        GameObject remotePlayerAvatar;
        remotePlayers.TryGetValue(player.username, out remotePlayerAvatar);
        Destroy(remotePlayerAvatar);
        remotePlayers.Remove(player.username);
        debugMessages.text = colyseusClient.GetClientList();
        updatePlayerList();
    }
                    
    public void OnPlayerChange(Player updatedPlayer)
    {
        // All player updates pass through here, including movement and interactions
        // though we need to handle those interactions differently. Keep this purely for movement!
        bool isLocal = updatedPlayer.username == localPlayer.username;
        
        Debug.Log("player id " + updatedPlayer.id + " is local: " + isLocal);

        GameObject remotePlayerAvatar;
        remotePlayers.TryGetValue(updatedPlayer.username, out remotePlayerAvatar);
        
        if (!isLocal)
        {
            if (remotePlayerAvatar != null)
            {
                // only move players with an avatar
                remotePlayerAvatar.GetComponent<RemotePlayerMovement>().NextPosition = new Vector3(
                    updatedPlayer.playerPosition.position.x, updatedPlayer.playerPosition.position.y,
                    updatedPlayer.playerPosition.position.z);
            }
        }
    }

    public void HandlePlayerInteraction(Player updatedPlayer, string interactionType)
    {
        // For all non-movement updates
        bool isLocal = updatedPlayer.username == localPlayer.username;
        if (!isLocal)
        {
            switch (interactionType)
            {
                case "interaction":
                    // show indicator
                    Debug.Log("Interaction update: " + updatedPlayer.interactionTarget.position.x + "," +
                              updatedPlayer.interactionTarget.position.y + "," +
                              updatedPlayer.interactionTarget.position.z);
                    Vector3 pos = Utils.NetworkPosToPosition(updatedPlayer.interactionTarget.position);
                    Quaternion rot = Utils.NetworkRotToRotation(updatedPlayer.interactionTarget.rotation);
                    ShowInteraction(pos, rot,
                        SimulationManager.GetInstance().GetColorForUsername(updatedPlayer.username), false);
                    break;
                case "celestialinteraction":
                    Debug.Log("remote player selected star");
                    // highlight star/ constellation
                    // TODO: Adjust how we create stars to make it possible to find the star from the network interaction
                    // this could be a simple rename, but need to check how constellation grouping works. Ideally we'll
                    // maintain a dict of stars by ID for easier lookups. 
                    SimulationManager.GetInstance().DataControllerComponent.GetStarById(updatedPlayer.celestialObjectTarget.uniqueId).HandleSelectStar();
                    break;
                case "locationpin":
                    // add / move player pin
                    Debug.Log("remote player pinned a location");
                    break;
                default:
                    break;
            }
        }
    }

    public void ShowInteraction(Vector3 pos, Quaternion rot, Color playerColor, bool isLocal)
    {
        if (interactionIndicator)
        {
            GameObject indicatorObj = Instantiate(interactionIndicator);
            indicatorObj.transform.localRotation = rot;
            indicatorObj.transform.position = pos;
            Utils.SetObjectColor(indicatorObj, playerColor);
            // Experimental code for Lat/Lng
            float radius = 5;
            if (GameObject.Find("Earth"))
            {
                radius = GameObject.Find("Earth").transform.localScale.x / 2;
            }
            Debug.Log("lat lng " + Utils.LatLngFromPosition(pos, radius));
            StartCoroutine(selfDestruct(indicatorObj));
        }
        if (isLocal)
        {
            // now need to broadcast to remotes
            colyseusClient.SendNetworkTransformUpdate(pos, rot, "interaction");
            string interactionInfo = "local interaction P:" +
    pos.ToString() + " R:" + rot.ToString();
            CCLogger.Log(CCLogger.EVENT_ADD_INTERACTION, interactionInfo);
        }
    }
    public void ShowCelestialObjectInteraction(string coName, string coGroup, string uniqueId, bool isLocal)
    {
        if (isLocal)
        {
            // now need to broadcast to remotes
            NetworkCelestialObject co = new NetworkCelestialObject
            {
                name = coName, group = coGroup, uniqueId = uniqueId
            };
            colyseusClient.SendCelestialInteraction(co);
            string interactionInfo = "local celestial interaction CO:" +
                                     coName + ", " + coGroup + ", " + uniqueId;
            CCLogger.Log(CCLogger.EVENT_ADD_INTERACTION, interactionInfo);
        }
        else
        {
            string interactionInfo = "remote celestial interaction CO:" +
                                     coName + ", " + coGroup + ", " + uniqueId;
        }
    }
    public void HandleMovementUpdate(Vector3 pos, Quaternion rot, bool isLocal)
    {
        if (isLocal)
        {
            // broadcast!
            colyseusClient.SendNetworkTransformUpdate(pos, rot, "movement");

            // log!
            string movementInfo = "local player moved to P:" +
                pos.ToString() + " R:" + rot.ToString();
            CCLogger.Log(CCLogger.EVENT_PLAYER_MOVE, movementInfo);
        }
        else
        {
            // move remote player
        }
    }

    IEnumerator selfDestruct(GameObject indicatorObj)
    {
        yield return new WaitForSeconds(3.0f);
        if (indicatorObj)
        {
            Destroy(indicatorObj);
        }
    }
    // called when the game is terminated
    void OnDisable()
    {
        Debug.Log("OnDisable");
        SceneManager.sceneLoaded -= OnSceneLoaded;
    }
    
}