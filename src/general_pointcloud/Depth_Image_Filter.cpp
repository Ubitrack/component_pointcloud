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
 * @ingroup pcl_components
 * @file
 * Componenten to filter a depth image
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

static log4cpp::Category& logger(log4cpp::Category::getInstance("Ubitrack.PointCloud.DepthMaptoPointCloud"));

namespace Ubitrack { namespace Components {


/**
 * @ingroup pcl_components
 */
        class PCL_FilterImage
                : public Dataflow::TriggerComponent
        {
        public:
            /**
             * UTQL component constructor.
             *
             * @param sName Unique name of the component.
             * @param subgraph UTQL subgraph
             */
            PCL_FilterImage( const std::string& sName, boost::shared_ptr< Graph::UTQLSubgraph > pConfig )
                    : Dataflow::TriggerComponent( sName, pConfig )
                    , m_inPortImage( "Input", *this )
                    , m_outPort( "Output", *this )
            {
            }

            /** Method that computes the result. */
            void compute( Measurement::Timestamp t )
            {

                Measurement::ImageMeasurement imageMea = m_inPortImage.get();
                Vision::Image::ImageFormatProperties fmt;
                imageMea->getFormatProperties(fmt);
                cv::Mat depthImage = imageMea->Mat();
                cv::Mat result = depthImage.clone();

                size_t count = depthImage.total();

                if(m_buffer.size() == 0) {
                    m_buffer = std::vector< std::vector<uint16_t> >(count , std::vector<uint16_t>() );
                }

                uint16_t* valueP = (uint16_t*)( depthImage.data );
                uint16_t* resultP = (uint16_t*)( result.data );
                for(int i=0;i<count;++i) {
                    m_buffer[i].push_back(valueP[i]);
                    std::sort(m_buffer[i].begin(), m_buffer[i].end());
                    if(m_buffer[i].size() % 2 == 0) {
                        // even
                        uint16_t p1 = m_buffer[i][m_buffer[i].size()/2-1];
                        uint16_t p2 = m_buffer[i][m_buffer[i].size()/2];
                        resultP[i] = (p1+p2)/2;
                    }
                    else {
                        // odd
                        resultP[i] = m_buffer[i][m_buffer[i].size()/2];
                    }


                }




                boost::shared_ptr< Vision::Image > pImage;
                pImage.reset(new Vision::Image(result, fmt));


                m_outPort.send( Measurement::ImageMeasurement( t, pImage ) );
            }

        protected:
            /** Input port A of the component. */
            Dataflow::TriggerInPort< Measurement::ImageMeasurement > m_inPortImage;


            /** Output port of the component. */
            Dataflow::TriggerOutPort< Measurement::ImageMeasurement > m_outPort;

            std::vector< std::vector<uint16_t> > m_buffer;



        };


        UBITRACK_REGISTER_COMPONENT( Dataflow::ComponentFactory* const cf ) {
        cf->registerComponent< PCL_FilterImage > ( "Depth_FilterImage" );


    }

} } // namespace Ubitrack::Components
