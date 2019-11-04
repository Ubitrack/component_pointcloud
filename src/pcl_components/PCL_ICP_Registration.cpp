/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */

/**
 * @ingroup driver_components
 * @file
 * ICP registration from pcl example
 *
 * 
 */

#include <utDataflow/TriggerComponent.h>
#include <utDataflow/TriggerInPort.h>
#include <utDataflow/TriggerOutPort.h>
#include <utDataflow/ComponentFactory.h>
#include <utMeasurement/Measurement.h>
#include <utVision/Image.h>
#include <utUtil/Logging.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

static log4cpp::Category& logger(log4cpp::Category::getInstance("Ubitrack.PointCloud.PCL_ICP"));

namespace Ubitrack { namespace Drivers {

        class PCL_ICP
                : public Dataflow::TriggerComponent {
        public:
            /**
             * UTQL component constructor.
             *
             * @param sName Unique name of the component.
             * @param subgraph UTQL subgraph
             */
            PCL_ICP(const std::string &sName, boost::shared_ptr<Graph::UTQLSubgraph> pConfig)
                    : Dataflow::TriggerComponent(sName, pConfig),
                      m_inPortPointCloud1("PC1", *this),
                      m_inPortPointCloud2("PC2", *this),
                      m_inPortPC1toPC2("PC1toPC2", *this),
                      m_outPort("Output", *this),
                      m_maxCorrespondenceDistance(0.05),
                      m_maximumIterations(50),
                      m_transformationEpsilon(1e-8),
                      m_euclideanFitnessEpsilon(1)
                      {

                        pConfig->m_DataflowAttributes.getAttributeData( "maxCorrespondenceDistance", m_maxCorrespondenceDistance );
                        pConfig->m_DataflowAttributes.getAttributeData( "maximumIterations", m_maximumIterations );
                        pConfig->m_DataflowAttributes.getAttributeData( "transformationEpsilon", m_transformationEpsilon );
                        pConfig->m_DataflowAttributes.getAttributeData( "euclideanFitnessEpsilon", m_euclideanFitnessEpsilon );
            }

            /** Method that computes the result. */
            void compute(Measurement::Timestamp t) {
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1 (new pcl::PointCloud<pcl::PointXYZ>);
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2 (new pcl::PointCloud<pcl::PointXYZ>);

                std::vector < Math::Vector< double, 3 > > pc1 = *m_inPortPointCloud1.get();
                std::vector < Math::Vector< double, 3 > > pc2 = *m_inPortPointCloud2.get();

                // Fill in the cloud_1 data
                cloud_1->width    = pc1.size();
                cloud_1->height   = 1;
                cloud_1->is_dense = false;
                cloud_1->points.resize (cloud_1->width * cloud_1->height);
                for (std::size_t i = 0; i < pc1.size(); ++i)
                {
                    cloud_1->points[i].x = static_cast<float>(pc1[i](0));
                    cloud_1->points[i].y = static_cast<float>(pc1[i](1));
                    cloud_1->points[i].z = static_cast<float>(pc1[i](2));
                }
                LOG4CPP_INFO(logger, "Point in cloud 1: " << pc1.size())

                // Fill in the cloud_2 data
                cloud_2->width    = pc2.size();
                cloud_2->height   = 1;
                cloud_2->is_dense = false;
                cloud_2->points.resize (cloud_2->width * cloud_2->height);
                for (std::size_t i = 0; i < pc2.size(); ++i)
                {
                    cloud_2->points[i].x = static_cast<float>(pc2[i](0));
                    cloud_2->points[i].y = static_cast<float>(pc2[i](1));
                    cloud_2->points[i].z = static_cast<float>(pc2[i](2));
                }
                LOG4CPP_INFO(logger, "Point in cloud 2: " << pc2.size())


                Math::Pose utGuess = *m_inPortPC1toPC2.get();
                Math::Matrix4x4d guessMatrix(utGuess);


                pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::Matrix4 guess;
                for(int i=0;i<16;++i){
                    guess(i) = static_cast<float>(guessMatrix(i));
                }


                pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;



                icp.setInputSource(cloud_1);
                icp.setInputTarget(cloud_2);

                // Set the max correspondence distance to 5cm (e.g., correspondences with higher distances will be ignored)
                icp.setMaxCorrespondenceDistance (m_maxCorrespondenceDistance);
                // Set the maximum number of iterations (criterion 1)
                icp.setMaximumIterations (m_maximumIterations);
                // Set the transformation epsilon (criterion 2)
                icp.setTransformationEpsilon (m_transformationEpsilon);
                // Set the euclidean distance difference epsilon (criterion 3)
                icp.setEuclideanFitnessEpsilon (m_euclideanFitnessEpsilon);

                pcl::PointCloud<pcl::PointXYZ> Final;
                icp.align(Final,guess);

                pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::Matrix4 finalPCLPose = icp.getFinalTransformation();
                Math::Matrix4x4d finalMatrix;

                for(int i=0;i<16;++i){
                    finalMatrix(i) = static_cast<double>(finalPCLPose(i));
                }
                Math::Pose finalUTPose(finalMatrix);

                LOG4CPP_INFO(logger, "ICP has converged:" << icp.hasConverged() << " score: " << icp.getFitnessScore());
                LOG4CPP_INFO(logger, "ICP initial guess: " << utGuess );
                LOG4CPP_INFO(logger, "ICP final transform: " << finalUTPose);




                m_outPort.send(Measurement::Pose(t, finalUTPose));
            }

        protected:
            /** Input port A of the component. */
            Dataflow::TriggerInPort<Measurement::PositionList> m_inPortPointCloud1;
            Dataflow::TriggerInPort<Measurement::PositionList> m_inPortPointCloud2;
            Dataflow::TriggerInPort<Measurement::Pose> m_inPortPC1toPC2;


            /** Output port of the component. */
            Dataflow::TriggerOutPort<Measurement::Pose> m_outPort;

            double m_maxCorrespondenceDistance;
            int m_maximumIterations;
            double m_transformationEpsilon;
            double m_euclideanFitnessEpsilon;


        };





        UBITRACK_REGISTER_COMPONENT( Dataflow::ComponentFactory* const cf ) {
            cf->registerComponent< PCL_ICP > ( "PCL_ICP" );


        }

} } // namespace Ubitrack::Drivers



