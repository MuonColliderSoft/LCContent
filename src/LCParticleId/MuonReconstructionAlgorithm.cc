/**
 *  @file   LCContent/src/LCParticleId/MuonReconstructionAlgorithm.cc
 * 
 *  @brief  Implementation of the muon reconstruction algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCParticleId/MuonReconstructionAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lc_content
{

MuonReconstructionAlgorithm::MuonReconstructionAlgorithm() :
    m_shouldClusterIsolatedHits(false),
    m_maxClusterCaloHits(30),
    m_minClusterOccupiedLayers(8),
    m_minClusterLayerSpan(8),
    m_nClusterLayersToFit(100),
    m_maxClusterFitChi2(20.f),
    m_maxDistanceToTrack(200.f),
    m_minTrackCandidateEnergy(7.f),
    m_minHelixClusterCosAngle(0.98f),
    m_nExpectedTracksPerCluster(1),
    m_nExpectedParentTracks(1),
    m_minHelixCaloHitCosAngle(0.95f),
    m_region1GenericDistance(3.f),
    m_region2GenericDistance(6.f),
    m_isolatedMinRegion1Hits(1),
    m_isolatedMaxRegion2Hits(0),
    m_maxGenericDistance(6.f),
    m_isolatedMaxGenericDistance(3.f),
    m_replaceCurrentClusterList(false),
    m_replaceCurrentPfoList(false),

    //LP: additional params
    m_insideOutSearch(false),
    m_forwardAngleDeg(25.f),
    m_transitionAngleDeg(35.f),
    m_minClusterOccupiedLayersFwd(5),
    m_minClusterLayerSpanFwd(5),
    m_maxHitsOverLayers(3.f),
    m_maxHitsOverLayersFwd(1.9f),
    m_maxDistanceToTrackFwd(200.f),
    m_maxDistanceToTrackTrn(200.f),
    m_maxRedTrackClusterChi2(50.f),
    m_maxRedTrackClusterChi2Fwd(50.f),
    m_maxRedTrackClusterChi2Trn(50.f),
    m_maxRedTrkCluChi2Direct(4.f),

    m_associateCaloHits(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::Run()
{
    std::string muonClusterListName;
    const ClusterList *pMuonClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_inputMuonCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::RunClusteringAlgorithm(*this, m_muonClusteringAlgName, pMuonClusterList, muonClusterListName));

    if (!pMuonClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AssociateMuonTracks(pMuonClusterList));
        if (m_associateCaloHits)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AddCaloHits(pMuonClusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateMuonPfos(pMuonClusterList));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TidyLists());

    return STATUS_CODE_SUCCESS;
}

//LP: muon propagation modified for compatibility with MUSIC detector -----------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::AssociateMuonTracks(const ClusterList *const pMuonClusterList) const
{
    const GeometryManager *const pGeometryManager(PandoraContentApi::GetGeometry(*this));
    const BFieldPlugin *const pBFieldPlugin(PandoraContentApi::GetPlugins(*this)->GetBFieldPlugin());
    const float innerBField(pBFieldPlugin->GetBField(CartesianVector(0.f, 0.f, 0.f)));

    const float coilMaxZ(std::fabs(pGeometryManager->GetSubDetector(COIL).GetOuterZCoordinate()));
    const float coilMidPointR(
        0.5f * (pGeometryManager->GetSubDetector(COIL).GetInnerRCoordinate() + pGeometryManager->GetSubDetector(COIL).GetOuterRCoordinate()));
    const float muonBarrelInnerR(std::fabs(pGeometryManager->GetSubDetector(MUON_BARREL).GetInnerRCoordinate()));
    const float muonBarrelBField(
        pBFieldPlugin->GetBField(CartesianVector(pGeometryManager->GetSubDetector(MUON_BARREL).GetInnerRCoordinate(), 0.f, 0.f)));
    const float muonEndCapBField(
        pBFieldPlugin->GetBField(CartesianVector(0.f, 0.f, std::fabs(pGeometryManager->GetSubDetector(MUON_ENDCAP).GetInnerZCoordinate()))));

    const TrackList *pTrackList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputTrackListName, pTrackList));

    if (m_insideOutSearch)
    {
        // Loop over all non-associated tracks in the current track list to find bestTrack
        for (TrackList::const_iterator iterT = pTrackList->begin(), iterTEnd = pTrackList->end(); iterT != iterTEnd; ++iterT)
        {
            const Track *const pTrack = *iterT;

            // Simple cuts on track properties
            if (pTrack->HasAssociatedCluster() || !pTrack->CanFormPfo())
                continue;

            if (!pTrack->GetDaughterList().empty())
                continue;

            if (pTrack->GetEnergyAtDca() < m_minTrackCandidateEnergy)
                continue;

            const bool isPositiveZ(pTrack->GetTrackStateAtCalorimeter().GetPosition().GetZ() > 0.f);

            // if (pTrack->IsProjectedToEndCap() && (pTrack->GetTrackStateAtCalorimeter().GetPosition().GetZ() * clusterInnerCentroid.GetZ() < 0.f))
            //     continue;

            // Extract track helix fit
            const Helix helix(pTrack->GetTrackStateAtCalorimeter().GetPosition(), pTrack->GetTrackStateAtCalorimeter().GetMomentum(),
                pTrack->GetCharge(), innerBField);

            // Propagate the track to the coil max Z
            CartesianVector muonFieldExitPoint(0.f, 0.f, 0.f);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMuonCoilExitPoint(helix, isPositiveZ, muonFieldExitPoint));
            CartesianVector muonFieldExitMomentum(helix.GetExtrapolatedMomentum(muonFieldExitPoint));
            const float helixCharge(helix.GetCharge());

            // Check if track is projected to muon endcap, barrel or transition
            bool isInEndcap(false), isInTransition(false);
            const float muonFieldExitR(
                std::sqrt(muonFieldExitPoint.GetX() * muonFieldExitPoint.GetX() + muonFieldExitPoint.GetY() * muonFieldExitPoint.GetY()));

            if (muonFieldExitR < coilMidPointR)
                isInEndcap = true; // the muon exits the solenoid towards the endcap
            else                   // the muon exits the solenoid from the barrel
            {
                const Helix outerCoilHelix(muonFieldExitPoint, muonFieldExitMomentum,
                    (muonBarrelBField < 0.f) ? -helix.GetCharge() : helix.GetCharge(), std::fabs(muonBarrelBField));

                CartesianVector muonOuterFieldExitPoint(0.f, 0.f, 0.f);
                float genericTime(std::numeric_limits<float>::max());

                const StatusCode statusCode(outerCoilHelix.GetPointInZ(
                    isPositiveZ ? coilMaxZ : -coilMaxZ, outerCoilHelix.GetReferencePoint(), muonOuterFieldExitPoint, genericTime));

                const float muonOuterFieldExitR(std::sqrt(muonOuterFieldExitPoint.GetX() * muonOuterFieldExitPoint.GetX() +
                                                          muonOuterFieldExitPoint.GetY() * muonOuterFieldExitPoint.GetY()));

                // the muon exits the solenoid from the barrel, then exits the outer field towards the endcap
                if (muonOuterFieldExitR < muonBarrelInnerR && STATUS_CODE_SUCCESS == statusCode)
                {
                    isInTransition = true;
                    muonFieldExitMomentum = outerCoilHelix.GetExtrapolatedMomentum(muonOuterFieldExitPoint);
                    muonFieldExitPoint = muonOuterFieldExitPoint;
                }
                // here it defaults to barrel = the muon exits the solenoid from the barrel and never exits the outer field
            }

            // Create helix that can be propagated in muon system, outside central detector
            float externalBField(muonBarrelBField);
            if (isInEndcap || isInTransition)
                externalBField = muonEndCapBField;

            // if(isInEndcap) std::cout << "Is endcap" << std::endl;
            // else if(isInTransition) std::cout << "Is transition" << std::endl;
            // else std::cout << "Is barrel" << std::endl;

            const Helix externalHelix(
                muonFieldExitPoint, muonFieldExitMomentum, (externalBField < 0.f) ? -helixCharge : helixCharge, std::fabs(externalBField));

            CartesianVector muonEntryPoint(0.f, 0.f, 0.f);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMuonEntryPoint(externalHelix, isPositiveZ, muonEntryPoint));

            const CartesianVector trackDir(pTrack->GetTrackStateAtStart().GetMomentum());
            float r, theta, phi;
            trackDir.GetSphericalCoordinates(r, phi, theta);
            theta = theta < 0.5 * 3.14159 ? theta : 3.14159 - theta;

            // float bestTrackEnergy(0.f);
            float bestDistanceToTrack(m_maxDistanceToTrack);
            if (theta < m_forwardAngleDeg * 3.14159 / 180.)
                bestDistanceToTrack = m_maxDistanceToTrackFwd;
            else if (theta < m_transitionAngleDeg * 3.14159 / 180.)
                bestDistanceToTrack = m_maxDistanceToTrackTrn;

            const Cluster *pBestCluster(NULL);

            for (ClusterList::const_iterator iter = pMuonClusterList->begin(), iterEnd = pMuonClusterList->end(); iter != iterEnd; ++iter)
            {
                const Cluster *const pCluster = *iter;

                unsigned int minClusterOccupiedLayers(m_minClusterOccupiedLayers);
                unsigned int minClusterLayerSpan(m_minClusterLayerSpan);
                float maxHitsOverLayers(m_maxHitsOverLayers);

                const CartesianVector clusterInnerCentroid(this->getCentroidInMuonYoke(pCluster)); //modified due to inner centroid miscalculation
                clusterInnerCentroid.GetSphericalCoordinates(r, phi, theta);
                theta = theta < 0.5 * 3.14159 ? theta : 3.14159 - theta;

                if (theta < m_forwardAngleDeg * 3.14159 / 180.)
                {
                    minClusterOccupiedLayers = m_minClusterOccupiedLayersFwd;
                    minClusterLayerSpan = m_minClusterLayerSpanFwd;
                    maxHitsOverLayers = m_maxHitsOverLayersFwd;
                }

                // Simple cuts on cluster properties
                if (pCluster->GetNCaloHits() > m_maxClusterCaloHits)
                    continue;

                if (pCluster->GetOrderedCaloHitList().size() < minClusterOccupiedLayers)
                    continue;

                if ((pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() + 1) < minClusterLayerSpan)
                    continue;

                if (pCluster->GetNCaloHits() / (pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() + 1) > maxHitsOverLayers)
                    continue;

                // Get direction of the cluster
                ClusterFitResult clusterFitResult;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitStart(pCluster, m_nClusterLayersToFit, clusterFitResult));

                if (!clusterFitResult.IsFitSuccessful())
                    continue;

                if (clusterFitResult.GetChi2() > m_maxClusterFitChi2)
                    continue;

                //const CartesianVector clusterInnerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
                // const bool isPositiveZ(clusterInnerCentroid.GetZ() > 0.f);

                const CartesianVector helixDirection(externalHelix.GetExtrapolatedMomentum(muonEntryPoint).GetUnitVector());
                const float helixClusterCosAngle(helixDirection.GetCosOpeningAngle(clusterFitResult.GetDirection()));

                if (helixClusterCosAngle < m_minHelixClusterCosAngle)
                    continue;

                float reducedTrackClusterChi2(this->GetChi2OfTrackClusterFit(pCluster, externalHelix));
                if (reducedTrackClusterChi2 > m_maxRedTrackClusterChi2)
                    continue;

                // Calculate separation of helix and cluster inner centroid
                CartesianVector helixSeparation(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, externalHelix.GetDistanceToPoint(clusterInnerCentroid, helixSeparation));

                const float distanceToTrack(helixSeparation.GetZ());

                if (distanceToTrack < bestDistanceToTrack)
                {
                    pBestCluster = pCluster;
                    bestDistanceToTrack = distanceToTrack;
                }
            }

            if (NULL != pBestCluster)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddTrackClusterAssociation(*this, pTrack, pBestCluster));
            }
        }
    }
    else
    {
        for (ClusterList::const_iterator iter = pMuonClusterList->begin(), iterEnd = pMuonClusterList->end(); iter != iterEnd; ++iter)
        {
            const Cluster *const pCluster = *iter;

            unsigned int minClusterOccupiedLayers(m_minClusterOccupiedLayers);
            unsigned int minClusterLayerSpan(m_minClusterLayerSpan);
            float maxHitsOverLayers(m_maxHitsOverLayers);
            float r, theta, phi;
            const CartesianVector clusterInnerCentroid(this->getCentroidInMuonYoke(pCluster)); //modified due to inner centroid miscalculation
            clusterInnerCentroid.GetSphericalCoordinates(r, phi, theta);
            theta = theta < 0.5 * 3.14159 ? theta : 3.14159 - theta;

            if (theta < m_forwardAngleDeg * 3.14159 / 180.)
            {
                minClusterOccupiedLayers = m_minClusterOccupiedLayersFwd;
                minClusterLayerSpan = m_minClusterLayerSpanFwd;
                maxHitsOverLayers = m_maxHitsOverLayersFwd;
            }

            // Simple cuts on cluster properties
            if (pCluster->GetNCaloHits() > m_maxClusterCaloHits)
                continue;

            if (pCluster->GetOrderedCaloHitList().size() < minClusterOccupiedLayers)
                continue;

            if ((pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() + 1) < minClusterLayerSpan)
                continue;

            if (pCluster->GetNCaloHits() / (pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() + 1) > maxHitsOverLayers)
                continue;

            // Get direction of the cluster
            ClusterFitResult clusterFitResult;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitStart(pCluster, m_nClusterLayersToFit, clusterFitResult));

            if (!clusterFitResult.IsFitSuccessful())
                continue;

            //const CartesianVector clusterInnerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
            const bool isPositiveZ(clusterInnerCentroid.GetZ() > 0.f);

            // Loop over all non-associated tracks in the current track list to find bestTrack
            const Track *pBestTrack(NULL);
            float bestTrackEnergy(0.f);
            //float bestDistanceToTrack(1000.f);
            float bestTrackClusterChi2(1000.f);

            for (TrackList::const_iterator iterT = pTrackList->begin(), iterTEnd = pTrackList->end(); iterT != iterTEnd; ++iterT)
            {
                const Track *const pTrack = *iterT;

                // Simple cuts on track properties
                if (pTrack->HasAssociatedCluster() || !pTrack->CanFormPfo())
                    continue;

                if (!pTrack->GetDaughterList().empty())
                    continue;

                if (pTrack->GetEnergyAtDca() < m_minTrackCandidateEnergy)
                    continue;

                if (pTrack->IsProjectedToEndCap() && (pTrack->GetTrackStateAtCalorimeter().GetPosition().GetZ() * clusterInnerCentroid.GetZ() < 0.f))
                    continue;

                // Extract track helix fit
                const Helix helix(pTrack->GetTrackStateAtCalorimeter().GetPosition(), pTrack->GetTrackStateAtCalorimeter().GetMomentum(),
                    pTrack->GetCharge(), innerBField);

                // Propagate the track to the coil max Z
                CartesianVector muonFieldExitPoint(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMuonCoilExitPoint(helix, isPositiveZ, muonFieldExitPoint));
                CartesianVector muonFieldExitMomentum(helix.GetExtrapolatedMomentum(muonFieldExitPoint));
                const float helixCharge(helix.GetCharge());

                // Check if track is projected to muon endcap, barrel or transition
                bool isInEndcap(false), isInTransition(false);
                const float muonFieldExitR(
                    std::sqrt(muonFieldExitPoint.GetX() * muonFieldExitPoint.GetX() + muonFieldExitPoint.GetY() * muonFieldExitPoint.GetY()));

                if (muonFieldExitR < coilMidPointR)
                    isInEndcap = true; // the muon exits the solenoid towards the endcap
                else                   // the muon exits the solenoid from the barrel
                {
                    const Helix outerCoilHelix(muonFieldExitPoint, muonFieldExitMomentum,
                        (muonBarrelBField < 0.f) ? -helix.GetCharge() : helix.GetCharge(), std::fabs(muonBarrelBField));

                    CartesianVector muonOuterFieldExitPoint(0.f, 0.f, 0.f);
                    float genericTime(std::numeric_limits<float>::max());

                    const StatusCode statusCode(outerCoilHelix.GetPointInZ(
                        isPositiveZ ? coilMaxZ : -coilMaxZ, outerCoilHelix.GetReferencePoint(), muonOuterFieldExitPoint, genericTime));

                    const float muonOuterFieldExitR(std::sqrt(muonOuterFieldExitPoint.GetX() * muonOuterFieldExitPoint.GetX() +
                                                              muonOuterFieldExitPoint.GetY() * muonOuterFieldExitPoint.GetY()));

                    // the muon exits the solenoid from the barrel, then exits the outer field towards the endcap
                    if (muonOuterFieldExitR < muonBarrelInnerR && STATUS_CODE_SUCCESS == statusCode)
                    {
                        isInTransition = true;
                        muonFieldExitMomentum = outerCoilHelix.GetExtrapolatedMomentum(muonOuterFieldExitPoint);
                        muonFieldExitPoint = muonOuterFieldExitPoint;
                    }
                    // here it defaults to barrel = the muon exits the solenoid from the barrel and never exits the outer field
                }

                // Create helix that can be propagated in muon system, outside central detector
                float externalBField(muonBarrelBField);
                if (isInEndcap || isInTransition)
                    externalBField = muonEndCapBField;

                const Helix externalHelix(muonFieldExitPoint, muonFieldExitMomentum, (externalBField < 0.f) ? -helixCharge : helixCharge,
                    std::fabs(externalBField));

                CartesianVector muonEntryPoint(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMuonEntryPoint(externalHelix, isPositiveZ, muonEntryPoint));

                const CartesianVector helixDirection(externalHelix.GetExtrapolatedMomentum(muonEntryPoint).GetUnitVector());
                const float helixClusterCosAngle(helixDirection.GetCosOpeningAngle(clusterFitResult.GetDirection()));

                if (helixClusterCosAngle < m_minHelixClusterCosAngle)
                    continue;

                // Calculate separation of helix and cluster inner centroid
                CartesianVector helixSeparation(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, externalHelix.GetDistanceToPoint(clusterInnerCentroid, helixSeparation));

                const float distanceToTrack(helixSeparation.GetZ());

                const CartesianVector trackDir(pTrack->GetTrackStateAtStart().GetMomentum());
                trackDir.GetSphericalCoordinates(r, phi, theta);
                theta = theta < 0.5 * 3.14159 ? theta : 3.14159 - theta;

                float reducedTrackClusterChi2(this->GetChi2OfTrackClusterFit(pCluster, externalHelix));

                bool passesCuts(this->PassTrackClusterCuts(clusterFitResult.GetChi2(), reducedTrackClusterChi2, distanceToTrack, theta));
                if (!passesCuts) continue;

                if ((std::fabs(reducedTrackClusterChi2-bestTrackClusterChi2)/bestTrackClusterChi2 < 0.1 && pTrack->GetEnergyAtDca() > bestTrackEnergy)
                    || reducedTrackClusterChi2 < 0.9 * bestTrackClusterChi2)
                {
                    pBestTrack = pTrack;
                    bestTrackClusterChi2 = reducedTrackClusterChi2;
                    bestTrackEnergy = pTrack->GetEnergyAtDca();
                    //std::cout << "Match! Theta = " << theta << std::endl;
                }
            }

            if (NULL != pBestTrack)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddTrackClusterAssociation(*this, pBestTrack, pCluster));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::GetMuonEntryPoint(const Helix &helix, const bool isPositiveZ, CartesianVector &muonEntryPoint) const
{
    const GeometryManager *const pGeometryManager(PandoraContentApi::GetGeometry(*this));
    const float muonEndCapInnerZ(std::fabs(pGeometryManager->GetSubDetector(MUON_ENDCAP).GetInnerZCoordinate()));

    float minGenericTime(std::numeric_limits<float>::max());
    const CartesianVector &referencePoint(helix.GetReferencePoint());

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        helix.GetPointInZ(isPositiveZ ? muonEndCapInnerZ : -muonEndCapInnerZ, referencePoint, muonEntryPoint, minGenericTime));

    const SubDetector &muonBarrel(pGeometryManager->GetSubDetector(MUON_BARREL));
    const unsigned int muonBarrelInnerSymmetry(muonBarrel.GetInnerSymmetryOrder());
    const float muonBarrelInnerPhi(muonBarrel.GetInnerPhiCoordinate());
    const float muonBarrelInnerR(muonBarrel.GetInnerRCoordinate());

    if (muonBarrelInnerSymmetry > 0)
    {
        const float pi(std::acos(-1.f));
        const float twopi_n = 2.f * pi / (static_cast<float>(muonBarrelInnerSymmetry));

        for (unsigned int i = 0; i < muonBarrelInnerSymmetry; ++i)
        {
            const float phi(twopi_n * static_cast<float>(i) + muonBarrelInnerPhi);

            CartesianVector barrelEntryPoint(0.f, 0.f, 0.f);
            float genericTime(std::numeric_limits<float>::max());

            const StatusCode statusCode(helix.GetPointInXY(muonBarrelInnerR * std::cos(phi), muonBarrelInnerR * std::sin(phi),
                std::cos(phi + 0.5f * pi), std::sin(phi + 0.5f * pi), referencePoint, barrelEntryPoint, genericTime));

            if ((STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
                minGenericTime = genericTime;
                muonEntryPoint = barrelEntryPoint;
            }
        }
    }
    else
    {
        CartesianVector barrelEntryPoint(0.f, 0.f, 0.f);
        float genericTime(std::numeric_limits<float>::max());

        const StatusCode statusCode(helix.GetPointOnCircle(muonBarrelInnerR, referencePoint, barrelEntryPoint, genericTime));

        if ((STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
        {
            minGenericTime = genericTime;
            muonEntryPoint = barrelEntryPoint;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//LP: added for compatibility with MUSIC detector ------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::GetMuonCoilExitPoint(const Helix &helix, const bool isPositiveZ, CartesianVector &muonExitPoint) const
{
    const GeometryManager *const pGeometryManager(PandoraContentApi::GetGeometry(*this));
    const float coilMaxZ(std::fabs(pGeometryManager->GetSubDetector(COIL).GetOuterZCoordinate()));

    float minGenericTime(std::numeric_limits<float>::max());
    const CartesianVector &referencePoint(helix.GetReferencePoint());

    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, helix.GetPointInZ(isPositiveZ ? coilMaxZ : -coilMaxZ, referencePoint, muonExitPoint, minGenericTime));

    const SubDetector &coilMagnet(pGeometryManager->GetSubDetector(COIL));
    const float coilMidPointR(0.5 * (coilMagnet.GetInnerRCoordinate() + coilMagnet.GetOuterRCoordinate()));

    CartesianVector coilBarrelIntersection(0.f, 0.f, 0.f);
    float genericTime(std::numeric_limits<float>::max());

    const StatusCode statusCode(helix.GetPointOnCircle(coilMidPointR, referencePoint, coilBarrelIntersection, genericTime));

    if ((STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
    {
        minGenericTime = genericTime;
        muonExitPoint = coilBarrelIntersection;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector MuonReconstructionAlgorithm::getCentroidInMuonYoke(const Cluster *const pCluster) const
{
    const GeometryManager *const pGeometryManager(PandoraContentApi::GetGeometry(*this));
    const float muonEndCapInnerZ(std::fabs(pGeometryManager->GetSubDetector(MUON_ENDCAP).GetInnerZCoordinate()));
    const float muonBarrelInnerR(pGeometryManager->GetSubDetector(MUON_BARREL).GetInnerRCoordinate());

    CartesianVector centroid(0.f, 0.f, 0.f);
    bool found = false;

    for (unsigned int i = pCluster->GetInnerPseudoLayer(); i < pCluster->GetOuterPseudoLayer() + 1; i++)
    {
        centroid = pCluster->GetCentroid(i);
        found = true;

        if (std::fabs(centroid.GetZ()) >= muonEndCapInnerZ)
            break;
        else if (std::sqrt(centroid.GetX() * centroid.GetX() + centroid.GetY() * centroid.GetY()) >= muonBarrelInnerR)
            break;
        else
            found = false;
    }

    if (found)
        return centroid;
    else
        return pCluster->GetCentroid(pCluster->GetInnerPseudoLayer());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MuonReconstructionAlgorithm::GetChi2OfTrackClusterFit(const Cluster *const pCluster, const Helix &helix) const
{
    float sumSquares(0.f);
    int nHits(0);
    CartesianVector helixD(0.f, 0.f, 0.f);

    const OrderedCaloHitList &orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    // loop over cluster hits, sum the normalized distances from the track helix
    for (const OrderedCaloHitList::value_type &layerIter : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *layerIter.second)
        {
            const StatusCode statusCode(helix.GetDistanceToPoint(pCaloHit->GetPositionVector(), helixD));
            if (STATUS_CODE_SUCCESS != statusCode)
                continue;
            float cellScale(pCaloHit->GetCellLengthScale());
            sumSquares += (helixD.GetZ() * helixD.GetZ()) / (cellScale * cellScale);
            nHits++;
        }
    }

    // return reduced chi2
    return sumSquares / float(nHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MuonReconstructionAlgorithm::PassTrackClusterCuts(const float &clusterFitChi2, const float &trackClusterChi2, const float &trackClusterD, float &trkTheta) const
{
    float maxDistanceToTrack(m_maxDistanceToTrack);
    float maxRedTrackClusterChi2(m_maxRedTrackClusterChi2);
    if (trkTheta < m_forwardAngleDeg * 3.14159 / 180.)
    {
        maxDistanceToTrack = m_maxDistanceToTrackFwd;
        maxRedTrackClusterChi2 = m_maxRedTrackClusterChi2Fwd;
    }
    else if (trkTheta < m_transitionAngleDeg * 3.14159 / 180.)
    {
        maxDistanceToTrack = m_maxDistanceToTrackTrn;
        maxRedTrackClusterChi2 = m_maxRedTrackClusterChi2Trn;
    }

    /*if (trackClusterChi2 < 300. && clusterFitChi2 < 300.)
    {
        std::cout << "Chi2 cluster = " << clusterFitChi2 << " / " << m_maxClusterFitChi2 << std::endl;
        std::cout << "Chi2 trk-clus = " << trackClusterChi2 << " / " << maxRedTrackClusterChi2 << std::endl;
        std::cout << "Th = " << trkTheta << "     Dist = " << trackClusterD << std::endl;
    }*/

    if (clusterFitChi2 > m_maxClusterFitChi2)
        return false;

    if (trackClusterChi2 < m_maxRedTrkCluChi2Direct)
        return true;

    if (trackClusterD < m_maxTrkCluDistanceDirect && trackClusterChi2 < m_maxRedTrackClusterChi2Close)
        return true;

    if (trackClusterChi2 < maxRedTrackClusterChi2 && trackClusterD < maxDistanceToTrack)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::AddCaloHits(const ClusterList *const pMuonClusterList) const
{
    const GeometryManager *const pGeometryManager(PandoraContentApi::GetGeometry(*this));
    const float hCalEndCapInnerR(pGeometryManager->GetSubDetector(HCAL_ENDCAP).GetInnerRCoordinate());
    const float eCalEndCapInnerR(pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetInnerRCoordinate());

    const SubDetector &coilMagnet(pGeometryManager->GetSubDetector(COIL));
    const float coilMidPointR(0.5 * (coilMagnet.GetInnerRCoordinate() + coilMagnet.GetOuterRCoordinate()));
    const float coilMaxZ(std::fabs(coilMagnet.GetOuterZCoordinate()));
    const float muonBarrelInnerR(std::fabs(pGeometryManager->GetSubDetector(MUON_BARREL).GetInnerRCoordinate()));

    const BFieldPlugin *const pBFieldPlugin(PandoraContentApi::GetPlugins(*this)->GetBFieldPlugin());
    const float muonBarrelBField(
        pBFieldPlugin->GetBField(CartesianVector(pGeometryManager->GetSubDetector(MUON_BARREL).GetInnerRCoordinate(), 0.f, 0.f)));
    const float muonEndCapBField(
        pBFieldPlugin->GetBField(CartesianVector(0.f, 0.f, std::fabs(pGeometryManager->GetSubDetector(MUON_ENDCAP).GetInnerZCoordinate()))));

    const float innerBField(pBFieldPlugin->GetBField(CartesianVector(0.f, 0.f, 0.f)));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    OrderedCaloHitList orderedCaloHitList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

    for (ClusterList::const_iterator clusterIter = pMuonClusterList->begin(), clusterIterEnd = pMuonClusterList->end();
        clusterIter != clusterIterEnd; ++clusterIter)
    {
        const Cluster *const pCluster = *clusterIter;

        // Check track associations
        const TrackList &trackList(pCluster->GetAssociatedTrackList());

        if (trackList.size() != m_nExpectedTracksPerCluster)
            continue;

        const Track *const pTrack = *(trackList.begin());
        //const Helix helix(pTrack->GetTrackStateAtCalorimeter().GetPosition(), pTrack->GetTrackStateAtCalorimeter().GetMomentum(), pTrack->GetCharge(), innerBField);
        const Helix *helix = new Helix(pTrack->GetTrackStateAtCalorimeter().GetPosition(),
            pTrack->GetTrackStateAtCalorimeter().GetMomentum(), pTrack->GetCharge(), innerBField);
        const bool isPositiveZ(pTrack->GetTrackStateAtCalorimeter().GetPosition().GetZ() > 0.f);

        // Propagate the track to the coil max Z
        CartesianVector muonFieldExitPoint(0.f, 0.f, 0.f);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMuonCoilExitPoint(*helix, isPositiveZ, muonFieldExitPoint));
        CartesianVector muonFieldExitMomentum(helix->GetExtrapolatedMomentum(muonFieldExitPoint));
        const float helixCharge(helix->GetCharge());
        // Check if track is projected to muon endcap, barrel or transition
        bool isInEndcap(false), isInTransition(false);
        const float muonFieldExitR(
            std::sqrt(muonFieldExitPoint.GetX() * muonFieldExitPoint.GetX() + muonFieldExitPoint.GetY() * muonFieldExitPoint.GetY()));

        const Helix *outerCoilHelix = nullptr;

        if (muonFieldExitR < coilMidPointR)
            isInEndcap = true; // the muon exits the solenoid towards the endcap
        else                   // the muon exits the solenoid from the barrel
        {
            outerCoilHelix = new Helix(muonFieldExitPoint, muonFieldExitMomentum,
                (muonBarrelBField < 0.f) ? -helix->GetCharge() : helix->GetCharge(), std::fabs(muonBarrelBField));

            CartesianVector muonOuterFieldExitPoint(0.f, 0.f, 0.f);
            float genericTime(std::numeric_limits<float>::max());

            const StatusCode statusCode(outerCoilHelix->GetPointInZ(
                isPositiveZ ? coilMaxZ : -coilMaxZ, outerCoilHelix->GetReferencePoint(), muonOuterFieldExitPoint, genericTime));

            const float muonOuterFieldExitR(std::sqrt(muonOuterFieldExitPoint.GetX() * muonOuterFieldExitPoint.GetX() +
                                                      muonOuterFieldExitPoint.GetY() * muonOuterFieldExitPoint.GetY()));

            // the muon exits the solenoid from the barrel, then exits the outer field towards the endcap
            if (muonOuterFieldExitR < muonBarrelInnerR && STATUS_CODE_SUCCESS == statusCode)
            {
                isInTransition = true;
                muonFieldExitMomentum = outerCoilHelix->GetExtrapolatedMomentum(muonOuterFieldExitPoint);
                muonFieldExitPoint = muonOuterFieldExitPoint;
            }
            // here it defaults to barrel = the muon exits the solenoid from the barrel and never exits the outer field
        }

        // Create helix that can be propagated in muon system, outside central detector
        float externalBField(muonBarrelBField);
        if (isInEndcap || isInTransition)
            externalBField = muonEndCapBField;

        // if(isInEndcap) std::cout << "Is endcap" << std::endl;
        // else if(isInTransition) std::cout << "Is transition" << std::endl;
        // else std::cout << "Is barrel" << std::endl;

        const Helix *externalHelix =
            new Helix(muonFieldExitPoint, muonFieldExitMomentum, (externalBField < 0.f) ? -helixCharge : helixCharge, std::fabs(externalBField));

        const Helix *refHelix = nullptr;

        for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
            layerIter != layerIterEnd; ++layerIter)
        {
            TrackDistanceInfoVector trackDistanceInfoVector;
            unsigned int nHitsInRegion1(0), nHitsInRegion2(0);

            for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end(); hitIter != hitIterEnd; ++hitIter)
            {
                const CaloHit *const pCaloHit = *hitIter;

                if ((!m_shouldClusterIsolatedHits && pCaloHit->IsIsolated()) || !PandoraContentApi::IsAvailable(*this, pCaloHit))
                    continue;

                const CartesianVector &caloHitPosition(pCaloHit->GetPositionVector());

                if ((caloHitPosition.GetZ() < 0 && isPositiveZ) || (caloHitPosition.GetZ() > 0 && !isPositiveZ))
                    continue;

                //choose helix segment based on calo hit position
                if (std::fabs(caloHitPosition.GetZ()) > coilMaxZ)
                {
                    if (isInTransition || isInEndcap)
                        refHelix = externalHelix;
                    else
                        continue;
                }
                else
                {
                    if (std::sqrt(std::pow(caloHitPosition.GetX(), 2) + std::pow(caloHitPosition.GetY(), 2)) > coilMidPointR)
                        refHelix = outerCoilHelix;
                    else
                        refHelix = helix;
                }

                if (!refHelix)
                    continue;

                const CartesianVector helixDirection(refHelix->GetExtrapolatedMomentum(caloHitPosition).GetUnitVector());

                if (pCaloHit->GetExpectedDirection().GetCosOpeningAngle(helixDirection) < m_minHelixCaloHitCosAngle)
                    continue;

                if (ENDCAP == pCaloHit->GetHitRegion())
                {
                    CartesianVector intersectionPoint(0.f, 0.f, 0.f);
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                        refHelix->GetPointInZ(caloHitPosition.GetZ(), refHelix->GetReferencePoint(), intersectionPoint));

                    const float helixR(
                        std::sqrt(intersectionPoint.GetX() * intersectionPoint.GetX() + intersectionPoint.GetY() * intersectionPoint.GetY()));

                    if ((HCAL == pCaloHit->GetHitType()) && (helixR < hCalEndCapInnerR))
                        continue;

                    if ((ECAL == pCaloHit->GetHitType()) && (helixR < eCalEndCapInnerR))
                        continue;
                }

                CartesianVector helixSeparation(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, refHelix->GetDistanceToPoint(caloHitPosition, helixSeparation));

                const float cellLengthScale(pCaloHit->GetCellLengthScale());

                if (cellLengthScale < std::numeric_limits<float>::epsilon())
                    continue;

                const float genericDistance(helixSeparation.GetZ() / cellLengthScale);
                trackDistanceInfoVector.push_back(TrackDistanceInfo(pCaloHit, genericDistance));

                if (genericDistance < m_region1GenericDistance)
                {
                    ++nHitsInRegion1;
                }
                else if (genericDistance < m_region2GenericDistance)
                {
                    ++nHitsInRegion2;
                }
            }

            const bool isIsolated((nHitsInRegion1 >= m_isolatedMinRegion1Hits) && (nHitsInRegion2 <= m_isolatedMaxRegion2Hits));
            std::sort(trackDistanceInfoVector.begin(), trackDistanceInfoVector.end(), MuonReconstructionAlgorithm::SortByDistanceToTrack);

            for (TrackDistanceInfoVector::const_iterator iter = trackDistanceInfoVector.begin(), iterEnd = trackDistanceInfoVector.end();
                iter != iterEnd; ++iter)
            {
                if ((iter->second > m_maxGenericDistance) || (isIsolated && (iter->second > m_isolatedMaxGenericDistance)))
                    break;

                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, iter->first));

                if (!isIsolated)
                    break;
            }
        }

        delete helix;
        delete externalHelix;
        delete outerCoilHelix;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::CreateMuonPfos(const ClusterList *const pMuonClusterList) const
{
    const PfoList *pPfoList = NULL;
    std::string pfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ClusterList::const_iterator iter = pMuonClusterList->begin(), iterEnd = pMuonClusterList->end(); iter != iterEnd; ++iter)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;

        const Cluster *const pCluster = *iter;
        pfoParameters.m_clusterList.push_back(pCluster);

        // Consider associated tracks
        const TrackList &trackList(pCluster->GetAssociatedTrackList());

        if (trackList.size() != m_nExpectedTracksPerCluster)
            continue;

        const Track *const pTrack = *(trackList.begin());
        pfoParameters.m_trackList.push_back(pTrack);

        // Examine track relationships
        const TrackList &parentTrackList(pTrack->GetParentList());

        if ((parentTrackList.size() > m_nExpectedParentTracks) || !pTrack->GetDaughterList().empty() || !pTrack->GetSiblingList().empty())
        {
            std::cout << "MuonReconstructionAlgorithm: invalid/unexpected track relationships for muon." << std::endl;
            continue;
        }

        if (!parentTrackList.empty())
        {
            pfoParameters.m_trackList.insert(pfoParameters.m_trackList.end(), parentTrackList.begin(), parentTrackList.end());
        }

        pfoParameters.m_charge = pTrack->GetCharge();
        pfoParameters.m_momentum = pTrack->GetMomentumAtDca();
        pfoParameters.m_particleId = (pfoParameters.m_charge.Get() > 0) ? MU_PLUS : MU_MINUS;
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy =
            std::sqrt(pfoParameters.m_mass.Get() * pfoParameters.m_mass.Get() + pfoParameters.m_momentum.Get().GetMagnitudeSquared());

        const ParticleFlowObject *pPfo(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pMuonClusterList->empty() && !pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputMuonClusterListName));

        if (m_replaceCurrentClusterList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputMuonClusterListName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputMuonPfoListName));

        if (m_replaceCurrentPfoList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputMuonPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::TidyLists() const
{
    // Make list of all tracks, clusters and calo hits in muon pfos
    TrackList pfoTrackList;
    CaloHitList pfoCaloHitList;
    ClusterList pfoClusterList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetPfoComponents(pfoTrackList, pfoCaloHitList, pfoClusterList));

    // Save the muon-removed track list
    const TrackList *pInputTrackList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputTrackListName, pInputTrackList));

    TrackList outputTrackList(*pInputTrackList);

    for (TrackList::const_iterator iter = pfoTrackList.begin(), iterEnd = pfoTrackList.end(); iter != iterEnd; ++iter)
    {
        TrackList::iterator outputIter = std::find(outputTrackList.begin(), outputTrackList.end(), *iter);

        if (outputTrackList.end() != outputIter)
            outputTrackList.erase(outputIter);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputTrackList, m_outputTrackListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Track>(*this, m_replacementTrackListName));

    // Save the muon-removed calo hit list
    const CaloHitList *pInputCaloHitList = NULL;
    const CaloHitList *pInputMuonCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pInputCaloHitList));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMuonCaloHitListName, pInputMuonCaloHitList));

    CaloHitList outputCaloHitList(*pInputCaloHitList);
    CaloHitList outputMuonCaloHitList(*pInputMuonCaloHitList);

    for (CaloHitList::const_iterator iter = pfoCaloHitList.begin(), iterEnd = pfoCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitList::iterator outputIter = std::find(outputCaloHitList.begin(), outputCaloHitList.end(), *iter);
        CaloHitList::iterator outputMuonIter = std::find(outputMuonCaloHitList.begin(), outputMuonCaloHitList.end(), *iter);

        if (outputCaloHitList.end() != outputIter)
            outputCaloHitList.erase(outputIter);

        if (outputMuonCaloHitList.end() != outputMuonIter)
            outputMuonCaloHitList.erase(outputMuonIter);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputCaloHitList, m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputMuonCaloHitList, m_outputMuonCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_replacementCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::GetPfoComponents(TrackList &pfoTrackList, CaloHitList &pfoCaloHitList, ClusterList &pfoClusterList) const
{
    const PfoList *pPfoList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_outputMuonPfoListName, pPfoList));

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo = *iter;
        const int particleId(pPfo->GetParticleId());

        if ((particleId != MU_MINUS) && (particleId != MU_PLUS))
            return STATUS_CODE_FAILURE;

        pfoTrackList.insert(pfoTrackList.end(), pPfo->GetTrackList().begin(), pPfo->GetTrackList().end());
        pfoClusterList.insert(pfoClusterList.end(), pPfo->GetClusterList().begin(), pPfo->GetClusterList().end());
    }

    for (ClusterList::const_iterator iter = pfoClusterList.begin(), iterEnd = pfoClusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(pfoCaloHitList);
        pfoCaloHitList.insert(pfoCaloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonReconstructionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Input lists
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputTrackListName", m_inputTrackListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputMuonCaloHitListName", m_inputMuonCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "MuonClusterFormation", m_muonClusteringAlgName));

    // Clustering
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldClusterIsolatedHits", m_shouldClusterIsolatedHits));

    // Cluster-track association
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterCaloHits", m_maxClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinClusterOccupiedLayers", m_minClusterOccupiedLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayerSpan", m_minClusterLayerSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NClusterLayersToFit", m_nClusterLayersToFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterFitChi2", m_maxClusterFitChi2));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToTrack", m_maxDistanceToTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinTrackCandidateEnergy", m_minTrackCandidateEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHelixClusterCosAngle", m_minHelixClusterCosAngle));

    // Addition of ecal/hcal hits
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NExpectedTracksPerCluster", m_nExpectedTracksPerCluster));

    if (0 == m_nExpectedTracksPerCluster)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NExpectedParentTracks", m_nExpectedParentTracks));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHelixCaloHitCosAngle", m_minHelixCaloHitCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Region1GenericDistance", m_region1GenericDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Region2GenericDistance", m_region2GenericDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IsolatedMinRegion1Hits", m_isolatedMinRegion1Hits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IsolatedMaxRegion2Hits", m_isolatedMaxRegion2Hits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGenericDistance", m_maxGenericDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IsolatedMaxGenericDistance", m_isolatedMaxGenericDistance));

    // Output lists
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTrackListName", m_outputTrackListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputMuonCaloHitListName", m_outputMuonCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputMuonClusterListName", m_outputMuonClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputMuonPfoListName", m_outputMuonPfoListName));

    // Current list management
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplacementTrackListName", m_replacementTrackListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplacementCaloHitListName", m_replacementCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentClusterList", m_replaceCurrentClusterList));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentPfoList", m_replaceCurrentPfoList));

    //LP: additional params
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InsideOutSearch", m_insideOutSearch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ForwardAngleDeg", m_forwardAngleDeg));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransitionAngleDeg", m_transitionAngleDeg));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinClusterOccupiedLayersFwd", m_minClusterOccupiedLayersFwd));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayerSpanFwd", m_minClusterLayerSpanFwd));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxHitsOverLayers", m_maxHitsOverLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxHitsOverLayersFwd", m_maxHitsOverLayersFwd));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToTrackFwd", m_maxDistanceToTrackFwd));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToTrackTrn", m_maxDistanceToTrackTrn));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AssociateCaloHits", m_associateCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxReducedTrackClusterChi2", m_maxRedTrackClusterChi2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxReducedTrackClusterChi2Fwd", m_maxRedTrackClusterChi2Fwd));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxReducedTrackClusterChi2Trn", m_maxRedTrackClusterChi2Trn));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxReducedTrackClusterChi2Direct", m_maxRedTrkCluChi2Direct));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTrackClusterDistanceDirect", m_maxTrkCluDistanceDirect));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxReducedTrackClusterChi2Close", m_maxRedTrackClusterChi2Close));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
